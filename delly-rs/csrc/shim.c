// Thin C shim exposing htslib macros / static inline functions to Rust.
#include <stdlib.h>
#include <string.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

// ---- SAM/BAM ---------------------------------------------------------------
int dl_sam_itr_next(htsFile *fp, hts_itr_t *itr, bam1_t *r) { return sam_itr_next(fp, itr, r); }
hts_itr_t *dl_sam_itr_queryi(const hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end) { return sam_itr_queryi(idx, tid, beg, end); }
int dl_hdr_n_targets(const sam_hdr_t *h) { return h->n_targets; }
const char *dl_hdr_target_name(const sam_hdr_t *h, int i) { return h->target_name[i]; }
uint32_t dl_hdr_target_len(const sam_hdr_t *h, int i) { return h->target_len[i]; }
const char *dl_hdr_text(const sam_hdr_t *h) { return h->text; }
size_t dl_hdr_l_text(const sam_hdr_t *h) { return h->l_text; }
int dl_bam_name2id(sam_hdr_t *h, const char *ref) { return bam_name2id(h, ref); }
int dl_hts_set_required_fields(htsFile *fp, int fields) { return hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, fields); }
int dl_hts_set_threads(htsFile *fp, int n) { return hts_set_threads(fp, n); }

// ---- VCF/BCF ---------------------------------------------------------------
int dl_bcf_hdr_nsamples(const bcf_hdr_t *h) { return bcf_hdr_nsamples(h); }
int dl_bcf_hdr_name2id(const bcf_hdr_t *h, const char *id) { return bcf_hdr_name2id(h, id); }
const char *dl_bcf_hdr_id2name(const bcf_hdr_t *h, int rid) { return bcf_hdr_id2name(h, rid); }
int dl_bcf_hdr_idinfo_exists_info(const bcf_hdr_t *h, int id) { return bcf_hdr_idinfo_exists(h, BCF_HL_INFO, id); }
int dl_bcf_hdr_id2int_id(const bcf_hdr_t *h, const char *id) { return bcf_hdr_id2int(h, BCF_DT_ID, id); }
void dl_bcf_set_rid(bcf1_t *r, int rid) { r->rid = rid; }
void dl_bcf_set_pos(bcf1_t *r, hts_pos_t pos) { r->pos = pos; }
void dl_bcf_set_qual(bcf1_t *r, float q) { r->qual = q; }
int dl_bcf_get_rid(const bcf1_t *r) { return r->rid; }
hts_pos_t dl_bcf_get_pos(const bcf1_t *r) { return r->pos; }
float dl_bcf_get_qual(const bcf1_t *r) { return r->qual; }
int dl_bcf_n_allele(const bcf1_t *r) { return r->n_allele; }
const char *dl_bcf_allele(const bcf1_t *r, int i) { return r->d.allele[i]; }
int dl_bcf_update_info_int32(const bcf_hdr_t *h, bcf1_t *r, const char *key, const int32_t *v, int n) { return bcf_update_info_int32(h, r, key, v, n); }
int dl_bcf_update_info_float(const bcf_hdr_t *h, bcf1_t *r, const char *key, const float *v, int n) { return bcf_update_info_float(h, r, key, v, n); }
int dl_bcf_update_info_flag(const bcf_hdr_t *h, bcf1_t *r, const char *key) { return bcf_update_info_flag(h, r, key, NULL, 1); }
int dl_bcf_update_info_string(const bcf_hdr_t *h, bcf1_t *r, const char *key, const char *s) { return bcf_update_info_string(h, r, key, s); }
int dl_bcf_update_format_int32(const bcf_hdr_t *h, bcf1_t *r, const char *key, const int32_t *v, int n) { return bcf_update_format_int32(h, r, key, v, n); }
int dl_bcf_update_genotypes(const bcf_hdr_t *h, bcf1_t *r, const int32_t *v, int n) { return bcf_update_genotypes(h, r, v, n); }
int dl_bcf_get_info_int32(const bcf_hdr_t *h, bcf1_t *r, const char *tag, int32_t **dst, int *ndst) { return bcf_get_info_int32(h, r, tag, dst, ndst); }
int dl_bcf_get_info_float(const bcf_hdr_t *h, bcf1_t *r, const char *tag, float **dst, int *ndst) { return bcf_get_info_float(h, r, tag, dst, ndst); }
int dl_bcf_get_info_string(const bcf_hdr_t *h, bcf1_t *r, const char *tag, char **dst, int *ndst) { return bcf_get_info_string(h, r, tag, dst, ndst); }
int dl_bcf_get_info_flag(const bcf_hdr_t *h, bcf1_t *r, const char *tag) { return bcf_get_info_flag(h, r, tag, 0, 0); }
int dl_bcf_unpack_info(bcf1_t *r) { return bcf_unpack(r, BCF_UN_INFO); }
