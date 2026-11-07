use std::env;

fn main() {
    // htslib location (Homebrew by default); must be the same htslib the reference delly links.
    let hts = env::var("HTSLIB_PREFIX").unwrap_or_else(|_| "/opt/homebrew".to_string());
    let inc = format!("{}/include", hts);
    let lib = format!("{}/lib", hts);

    cc::Build::new().file("csrc/shim.c").include(&inc).opt_level(3).compile("dellyshim");
    cc::Build::new()
        .cpp(true)
        .file("csrc/edlib.cpp")
        .file("csrc/edshim.cpp")
        .include("csrc")
        .flag("-std=c++17")
        .define("NDEBUG", None)
        .opt_level(3)
        .compile("edlib");

    println!("cargo:rustc-link-search=native={}", lib);
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib);
    println!("cargo:rustc-link-lib=dylib=hts");
    println!("cargo:rerun-if-changed=csrc");
    println!("cargo:rerun-if-env-changed=HTSLIB_PREFIX");
}
