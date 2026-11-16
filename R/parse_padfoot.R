# Padfoot SV annotation: parse padfoot/<combo>/annotated_svs.tsv and join it onto a caller's SV rows.
#
# Plain source()d helpers, no section registration. Sourced by templates/per_sample.qmd and
# tests/testthat/setup.R after R/parse_severus.R and before R/sections/*, so the SV section
# (R/sections/sv.R) calls into it.
#
# The format is taken from the pipeline's Padfoot module and its nf-test (no real output was on
# disk when this was written): 26 tab-separated columns,
#   SV ID, chr_1, pos_1, chr_2, pos_2, N_supp, Vaf, Strand, Score, Impact, Role in Cancer1,
#   Gene Name1, Pos_1, Strand_1, Repeat_Anno1, Role in Cancer2, Gene Name2, Pos_2, Strand_2,
#   Repeat_Anno2, Type, Microhomology, VNTR, Telomere_repeat, Aligned_pos(INS), Repeat_annot
# `pos_1` (a coordinate) and `Pos_1` (where in the gene: "intron3", "exon10", "downstream")
# differ only in case, which is why resolve_columns() is given exact names first.
#
# Padfoot's SV ID is not the same key for every caller: for Severus a breakend pair is reported
# under its `_1` record (`severus_BND1_1`), for SAVANA under the shared stem (`ID_1` for
# `ID_1_1`/`ID_1_2`). join_padfoot() therefore tries several keys and picks ONE for the whole
# table — whichever matches most rows — exactly as build_sv_table_from_vep() does.

suppressPackageStartupMessages({
  library(data.table)
})

# Columns join_padfoot() adds to an SV table, in display order
PADFOOT_COLS = c("padfoot_impact", "padfoot_score", "padfoot_role", "padfoot_genes",
                 "padfoot_type")

# Canonical -> candidate header names (exact names first; resolve_columns() then normalises)
PADFOOT_HEADER = list(
  sv_id   = c("SV ID", "SV_ID", "sv_id", "ID"),
  chr_1   = "chr_1",  pos_1 = "pos_1",
  chr_2   = "chr_2",  pos_2 = "pos_2",
  score   = "Score",  impact = "Impact", type = "Type",
  role_1  = c("Role in Cancer1", "Role in Cancer_1", "Role_in_Cancer1"),
  gene_1  = c("Gene Name1", "Gene Name_1", "Gene_Name1"),
  where_1 = "Pos_1",
  role_2  = c("Role in Cancer2", "Role in Cancer_2", "Role_in_Cancer2"),
  gene_2  = c("Gene Name2", "Gene Name_2", "Gene_Name2"),
  where_2 = "Pos_2"
)

# Cell values Padfoot (pandas) writes for "nothing here"
.padfoot_blank = function(x) {
  x = trimws(as.character(x))
  is.na(x) | x %in% c("", "NA", "NaN", "nan", "None", "none", ".", "-")
}

# Padfoot annotated_svs.tsv -> data.table(sv_id, padfoot_impact, padfoot_score, padfoot_role,
# padfoot_genes, padfoot_type, chrom_1, pos_1, chrom_2, pos_2); NULL when the file is absent,
# an empty table when it has no rows. The locus columns serve join_padfoot()'s fallback only.
parse_padfoot = function(f) {
  if (is.null(f) || length(f) == 0 || is.na(f) || !file.exists(f)) return(NULL)
  empty = data.table(sv_id = character(), padfoot_impact = character(),
                     padfoot_score = numeric(), padfoot_role = character(),
                     padfoot_genes = character(), padfoot_type = character(),
                     chrom_1 = character(), pos_1 = integer(),
                     chrom_2 = character(), pos_2 = integer())
  if (file.info(f)$size == 0) return(empty)

  # Everything as character, no quote handling: free-text annotation cells must not shift columns
  dt = tryCatch(fread(f, sep = "\t", quote = "", colClasses = "character",
                      check.names = FALSE),
                error = function(e) {
                  message("Padfoot: could not read ", f, ": ", conditionMessage(e))
                  NULL
                })
  if (is.null(dt) || nrow(dt) == 0) return(empty)

  cols = resolve_columns(dt, PADFOOT_HEADER)
  if (is.na(cols[["sv_id"]])) {
    message("Padfoot: no SV ID column in ", f, " - ignored.")
    return(empty)
  }
  get = function(k) {
    v = if (is.na(cols[[k]])) rep(NA_character_, nrow(dt)) else as.character(dt[[cols[[k]]]])
    v[.padfoot_blank(v)] = NA_character_
    trimws(v)
  }

  g1 = get("gene_1"); w1 = get("where_1"); g2 = get("gene_2"); w2 = get("where_2")
  r1 = get("role_1"); r2 = get("role_2")

  # "GENE (intron3)" per side; both sides, de-duplicated, so a deletion inside one gene reads once per position
  side_label = function(g, w) ifelse(is.na(g), NA_character_,
                                     ifelse(is.na(w), g, paste0(g, " (", w, ")")))
  l1 = side_label(g1, w1); l2 = side_label(g2, w2)
  genes = vapply(seq_len(nrow(dt)), function(i) {
    x = unique(c(l1[i], l2[i])); x = x[!is.na(x)]
    if (length(x) == 0) NA_character_ else paste(x, collapse = ", ")
  }, character(1))

  # Roles of both sides; a cell can itself hold "oncogene, fusion"
  roles = vapply(seq_len(nrow(dt)), function(i) {
    x = trimws(unlist(strsplit(c(r1[i], r2[i]), "[,;/]+")))
    x = unique(x[!is.na(x) & nzchar(x) & !.padfoot_blank(x)])
    if (length(x) == 0) NA_character_ else paste(x, collapse = ", ")
  }, character(1))

  chrom_or_na = function(x) { o = ensure_chr_prefix(x); o[is.na(x)] = NA_character_; o }

  data.table(
    sv_id          = get("sv_id"),
    padfoot_impact = get("impact"),
    padfoot_score  = suppressWarnings(as.numeric(get("score"))),
    padfoot_role   = roles,
    padfoot_genes  = genes,
    padfoot_type   = get("type"),
    chrom_1        = chrom_or_na(get("chr_1")),
    pos_1          = suppressWarnings(as.integer(get("pos_1"))),
    chrom_2        = chrom_or_na(get("chr_2")),
    pos_2          = suppressWarnings(as.integer(get("pos_2")))
  )
}

# Join Padfoot rows onto one caller's SV table (SV_TABLE_COLS; uses id, id_b when present, and
# the bedpe loci). Returns list(table, key, n_joined, n_rows):
#   table     sv with PADFOOT_COLS appended (NA where nothing joined)
#   key       "id" | "mate id" | "id stem" | "locus" | NA (nothing joined)
#   n_joined  rows that received a Padfoot row — a near-zero count is how a wrong key surfaces
# One key for the whole table (most matches wins; ties go id > mate id > id stem), so one side
# of the table cannot silently be joined differently from the other. Only when no id key
# matches at all does it fall back to the loci (both ends, either orientation, then the first
# end alone where that locus is unique in the Padfoot table).
join_padfoot = function(sv, pf) {
  out = copy(as.data.table(sv))
  out[, `:=`(padfoot_impact = NA_character_, padfoot_score = NA_real_,
             padfoot_role = NA_character_, padfoot_genes = NA_character_,
             padfoot_type = NA_character_)]
  res = list(table = out, key = NA_character_, n_joined = 0L, n_rows = nrow(out))
  if (is.null(pf) || nrow(pf) == 0 || nrow(out) == 0) return(res)

  pf = pf[!is.na(sv_id) & !duplicated(sv_id)]
  ids  = as.character(out$id)
  mate = if ("id_b" %in% names(out)) as.character(out$id_b) else rep(NA_character_, nrow(out))
  cand = list("id" = ids, "mate id" = mate, "id stem" = sub("_[12]$", "", ids))
  hits = vapply(cand, function(k) sum(!is.na(k) & k %in% pf$sv_id), integer(1))

  if (max(hits) > 0) {
    key = names(hits)[which.max(hits)]
    m   = match(cand[[key]], pf$sv_id)
  } else {
    key = "locus"
    lk  = function(c1, p1, c2, p2) paste(c1, p1, c2, p2, sep = ":")
    pf_key = lk(pf$chrom_1, pf$pos_1, pf$chrom_2, pf$pos_2)
    m = match(lk(out$chrom_a, out$pos_a, out$chrom_b, out$pos_b), pf_key)
    miss = is.na(m)
    m[miss] = match(lk(out$chrom_b, out$pos_b, out$chrom_a, out$pos_a)[miss], pf_key)
    # First end alone (an INS whose second end Padfoot leaves blank), only where unambiguous
    one = paste(pf$chrom_1, pf$pos_1, sep = ":")
    uniq_one = !(one %in% one[duplicated(one)])
    miss = is.na(m)
    m1 = match(paste(out$chrom_a, out$pos_a, sep = ":")[miss], one[uniq_one])
    m[miss] = which(uniq_one)[m1]
  }

  hit = !is.na(m)
  if (any(hit)) {
    for (nm in PADFOOT_COLS) set(out, which(hit), nm, pf[[nm]][m[hit]])
  }
  list(table = out, key = if (any(hit)) key else NA_character_,
       n_joined = sum(hit), n_rows = nrow(out))
}
