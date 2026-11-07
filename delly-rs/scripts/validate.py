#!/usr/bin/env python3
"""Compare callers on the bundled example and generated ordering/cap fixtures.

Usage: python3 scripts/validate.py <reference> [<rust-binary>] [--outdir DIR]
Requires samtools and bcftools. Artifacts/logs are retained for inspection.
The reference may also be a saved Rust binary when checking optimization parity.
"""
import argparse
import difflib
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT.parent / "delly/example"


def command(args, **kwargs):
    return subprocess.run([str(a) for a in args], check=True, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("rust", type=Path, nargs="?", default=ROOT / "target/release/delly")
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--comparison", choices=("bytes", "records"), default="bytes",
        help="bytes: compare BCF/CSI/dumps; records: compare sample order and every VCF record field")
    opts = parser.parse_args()
    out = opts.outdir or Path(tempfile.mkdtemp(prefix="delly-validate-"))
    out.mkdir(parents=True, exist_ok=True)
    out = out.resolve()
    print(f"Artifacts: {out}", flush=True)
    binaries = [("ref", opts.reference.resolve()), ("rs", opts.rust.resolve())]

    def call(binary, name, ref, bams, extra=(), threads=1, dump=False, batch=None, stdout=False):
        bcf = out / f"{name}.bcf"
        args = [binary, "sr", "-g", ref, "--sex", "none", "-h", str(threads)]
        if not stdout:
            args += ["-o", bcf]
        if dump:
            args += ["-d", out / f"{name}.dump.gz"]
        args += list(extra) + list(bams)
        env = dict(os.environ, DELLY_TIMING="1")
        env.pop("DELLY_BATCH_BASE", None)
        if batch is not None:
            env["DELLY_BATCH_BASE"] = str(batch)
        with (out / f"{name}.log").open("wb") as log:
            got = command(args, stdout=subprocess.PIPE, stderr=log, env=env)
        if opts.comparison == "records":
            path = "-" if stdout else str(bcf)
            source = {"input": got.stdout} if stdout else {}
            samples = subprocess.check_output(["bcftools", "query", "-l", path], **source)
            records = subprocess.check_output(["bcftools", "view", "--no-version", "-H", path], **source)
            (out / f"{name}.samples.txt").write_bytes(samples)
            (out / f"{name}.records.vcf").write_bytes(records)
            return samples, records
        if stdout:
            return (got.stdout,)
        values = [bcf.read_bytes(), Path(str(bcf) + ".csi").read_bytes()]
        if dump:
            values.append((out / f"{name}.dump.gz").read_bytes())
        return tuple(values)

    ref = EXAMPLE / "ref.fa"
    bam = EXAMPLE / "sr.bam"
    call(binaries[0][1], "seed", ref, [bam])
    vcf = subprocess.check_output(["bcftools", "view", "-Ov", str(out / "seed.bcf")], text=True)
    headers = [line for line in vcf.splitlines() if line.startswith("#")]
    record = next(line for line in vcf.splitlines() if not line.startswith("#"))
    tied = out / "tied.vcf"
    rows = []
    for i in range(256):
        fields = record.split("\t")
        fields[2] = f"DEL{i:08d}"
        rows.append("\t".join(fields))
    tied.write_text("\n".join(headers + rows) + "\n")

    # Two identical reference contigs with independent intra-chromosome reads,
    # plus actual inter-chromosome mate pairs copied from discordant DEL reads.
    sequence = "".join(line for line in ref.read_text().splitlines() if not line.startswith(">"))
    mixed_ref = out / "mixed.fa"
    mixed_ref.write_text(f">18\n{sequence}\n>chrB\n{sequence}\n")
    command(["samtools", "faidx", mixed_ref])
    header = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
    reads = subprocess.check_output(["samtools", "view", str(bam)], text=True).splitlines()
    sam = out / "mixed.sam"
    with sam.open("w") as stream:
        stream.write(header.replace("@SQ\tSN:18", f"@SQ\tSN:chrB\tLN:{len(sequence)}\n@SQ\tSN:18"))
        # chrB is tid 0, 18 is tid 1: BND primary coordinates use 18 (the later tid).
        for line in reads:
            stream.write(line + "\n")
            fields = line.split("\t")
            fields[0] = "mirror_" + fields[0]
            fields[2] = "chrB"
            if fields[6] == "18":
                fields[6] = "chrB"
            stream.write("\t".join(fields) + "\n")
        for line in reads:
            fields = line.split("\t")
            if abs(int(fields[8])) <= 5000 or int(fields[1]) & (0x100 | 0x800):
                continue
            fields[0] = "tra_" + fields[0]
            left = int(fields[3]) < int(fields[7])
            fields[2], fields[6] = ("18", "chrB") if left else ("chrB", "18")
            fields[8] = "0"
            stream.write("\t".join(fields) + "\n")
    mixed_bam = out / "mixed.bam"
    command(["samtools", "sort", "-o", mixed_bam, sam])
    command(["samtools", "index", mixed_bam])
    mixed_cram = out / "mixed.cram"
    command(["samtools", "view", "-C", "-T", mixed_ref, "-o", mixed_cram, mixed_bam])
    command(["samtools", "index", mixed_cram])
    mixed_sites = out / "mixed.vcf"
    mixed_headers = headers[:-1] + [f"##contig=<ID=chrB,length={len(sequence)}>", headers[-1]]
    mixed_rows = [record]
    for i in range(32):
        fields = record.split("\t")
        fields[2] = f"BND{i:08d}"
        info = dict(field.split("=", 1) if "=" in field else (field, None) for field in fields[7].split(";"))
        end = info["END"]
        fields[4] = f"{fields[3]}[chrB:{end}["
        info.update(SVTYPE="BND", END=fields[1], CHR2="chrB", POS2=end)
        info.pop("SVLEN", None)
        fields[7] = ";".join(key if value is None else f"{key}={value}" for key, value in info.items())
        mixed_rows.append("\t".join(fields))
    mixed_sites.write_text("\n".join(mixed_headers + mixed_rows) + "\n")

    cases = []
    for threads in (1, 4, 8):
        for dump in (False, True):
            cases.append((f"discovery-h{threads}-d{int(dump)}", ref, [bam], [], threads, dump, None, False))
    for cap in (0, 1, 5, 250):
        for threads in (1, 8):
            for dump in (False, True):
                cases.append((f"tied-a{cap}-h{threads}-d{int(dump)}", ref, [bam], ["-v", tied, "-a", str(cap)], threads, dump, 3, False))
    for cap in (1, 20, 250):
        for threads in (1, 8):
            for dump in (False, True):
                cases.append((f"bnd-a{cap}-h{threads}-d{int(dump)}", mixed_ref, [mixed_bam], ["-v", mixed_sites, "-a", str(cap)], threads, dump, 3, False))
    for dump in (False, True):
        cases.append((f"cram-d{int(dump)}", mixed_ref, [mixed_cram], ["-v", mixed_sites, "-a", "20"], 4, dump, 3, False))
    cases.append(("multisample", mixed_ref, [mixed_bam, mixed_bam], ["-v", mixed_sites, "-a", "5"], 4, True, 3, False))
    cases.append(("stdout", ref, [bam], [], 4, False, None, True))
    failures = []
    for name, genome, bams, extra, threads, dump, batch, stdout in cases:
        values = [call(binary, f"{name}.{impl}", genome, bams, extra, threads, dump, batch, stdout)
            for impl, binary in binaries]
        if values[0] != values[1]:
            failures.append(name)
            if opts.comparison == "records":
                for i, kind in enumerate(("samples", "records")):
                    diff = difflib.unified_diff(values[0][i].decode().splitlines(True),
                        values[1][i].decode().splitlines(True), fromfile="reference", tofile="rust")
                    (out / f"{name}.{kind}.diff").write_text("".join(diff))
            print(f"FAIL {name}", flush=True)
        else:
            count = f" ({len(values[0][1].splitlines())} rows)" if opts.comparison == "records" else ""
            print(f"PASS {name}{count}", flush=True)
    if failures:
        raise SystemExit(f"{len(failures)}/{len(cases)} cases differ: {', '.join(failures)}; inspect {out}")
    label = "identical VCF records and sample order" if opts.comparison == "records" else "byte-identical outputs"
    print(f"All {len(cases)} cases have {label}.")


if __name__ == "__main__":
    main()
