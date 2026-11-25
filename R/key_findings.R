# Clinical overview: key findings (driver summary) + genome-wide biomarkers.
#
# build_overview(variant_table, section_data, all_panels, ascat_pp, tmb_info, counts, ctx, hmf)
#   hmf: build_hmf()'s result (NULL without --hmf-*); the hmf_driver list's findings (kf_hmf_rows())
#   variant_table  FULL small-variant table (build_variant_table()), may be NULL
#   section_data   SECTION_DATA (sv, whatshap, cna, shatterseek, ecdna, reconplot, signatures)
#   all_panels     named list of gene panels (after the YAML round trip: unlist() fields)
#   ascat_pp       list(purity, ploidy) from parse_ascat_purityploidy()
#   tmb_info       compute_tmb() result: list(n_nonsyn, tmb, denominator_mb)
#   counts         list(small_variants = count_small_variants() result, n_sv, n_trans,
#                       n_intra_bnd, sv_caller)
#   ctx            the section ctx (see R/sections.R); reads $reference and $outputs$mode
# returns
#   list(findings    = data.table(gene, type, alteration, evidence, tier (integer), role,
#                                 panels, vaf, source_id,            <- the stable contract
#                                 genes,                             <- member symbols
#                                 dp, callers, locus, loh, in_panel, <- extra, display only
#                                 rescued, ev_tier),                 <- whitelist-rescued SVs
#                      A whitelist-rescued SV row has tier KF_TIER_RESCUED and `rescued` TRUE;
#                      `ev_tier` is the tier its evidence alone gives (= tier on other rows).
#                      `gene` is the row's group key: one gene, or "EWSR1–ATF1" for a
#                      rearrangement joining two listed genes, whose `genes` is
#                      "EWSR1,ATF1" (equal to `gene` on every other row). Rows are ordered
#                      so that a group's rows are contiguous: group best tier, group, tier.
#        biomarkers  = named list (purity, wgd, loh, tmb, small, sv, chromothripsis, ecdna,
#                      signatures) of list(id, label, value, detail, status ("ok"|"not_run"|
#                      "warn"), [items: character, signatures only]); never NULL, absent
#                      sections read value "not run",
#        per_panel_n = named integer (genes / groups with a tier 1-2 finding per panel;
#                      every panel named),
#        tumour_only = logical(1),
#        verdict     = parse_verdict_status() result or NULL (ClairS-TO germline tagging),
#        n_nonpanel  = list(kept, total))  # non-panel small-variant rows embedded / eligible
#
# What is a finding: an event in a gene of a user-chosen gene list ("panel") that carries
# evidence. There is no bundled pan-cancer driver list, by decision. Rows are computed for
# EVERY embedded panel at once and tagged with the panels whose gene set contains the gene
# (`panels`, comma-joined); assets/js/key_findings.js shows the rows of the ticked panels.
# A gene sits in a panel per table: small variants use the panel's SNV scope (genes_snv via
# panel_snv_genes()), SVs its SV scope (sv_panel_hit_table(), genes_sv), CNA the panel's
# gene-level CN rows (SECTION_DATA$cna$genes, built from every row with coordinates).
#
# Non-panel rows exist only so a browser-typed "Custom…" list has rows to match: small
# variants that are ClinVar (likely) pathogenic or HIGH impact, outside every panel, capped
# at KF_MAX_NONPANEL (strongest evidence first). `panels` is "" on those rows.
#
# Tiers (kf_tier()): 1 = strong evidence in a list gene, 2 = weaker or role-discordant,
# 3 = low confidence (SV merely near a gene, CN in a pseudoautosomal or IG/TR locus),
# "rescued" (KF_TIER_RESCUED) = an SV kept only by the Severus whitelist, hidden by default.

suppressPackageStartupMessages({
  library(data.table)
})

# Pseudoautosomal regions, 1-based inclusive, per reference (hg38: GRC; t2t: chm13v2.0_PAR.bed).
# A copy-number call here is unreliable: reads from the chrX and chrY copies map to one of them,
# so a PAR gene can read CN 0 with nothing deleted (CRLF2 on a real T2T sample did).
KF_PAR = list(
  hg38 = data.table(chrom = c("chrX", "chrX", "chrY", "chrY"),
                    start = c(10001, 155701383, 10001, 56887903),
                    end   = c(2781479, 156030895, 2781479, 57217415)),
  t2t  = data.table(chrom = c("chrX", "chrX", "chrY", "chrY"),
                    start = c(0, 153925834, 0, 62122809),
                    end   = c(2394410, 154259566, 2458320, 62460029))
)

# Immunoglobulin / T-cell receptor loci: physiologically rearranged (and deleted) in lymphocytes,
# so a copy-number loss there is not a driver call. SVs keep their tier: an IG partner is the point
KF_RECOMBINING_LOCI = "^(IG[HKL]|TR[ABDG])([^A-Z0-9]|$)|^TRA/D$"

KF_MAX_PROXIMITY_PER_GENE = 3L    # tier-3 "near GENE" SV rows kept per gene
KF_MAX_NONPANEL           = 500L  # non-panel small-variant rows embedded for Custom… lists
KF_FOCAL_SV_MAX           = 3e6   # a DEL/DUP up to this size counts when it covers a gene
# Two callers' calls are one event when both breakpoints agree within this distance (and
# within half the event's size, so that two neighbouring small deletions stay apart)
KF_SV_SAME_EVENT_BP       = 1000
# Tier of a whitelist-rescued SV finding (Severus INFO/WL_RESCUE): below tier 3, shown on request
KF_TIER_RESCUED           = 4L

# High-impact consequences that truncate / abolish the product (loss of function)
KF_LOF_TERMS = c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
                 "stop_gained", "frameshift_variant", "start_lost")

# Consequences that put a variant in a gene's coding sequence or splice sites: a ClinVar
# classification (variant-level) is credited to the gene rows carrying one of these
KF_CLINVAR_GENE_TERMS = c(KF_LOF_TERMS, "missense_variant", "synonymous_variant",
                          "stop_lost", "stop_retained_variant", "start_retained_variant",
                          "inframe_insertion", "inframe_deletion", "protein_altering_variant",
                          "coding_sequence_variant", "incomplete_terminal_codon_variant",
                          "splice_region_variant", "splice_donor_5th_base_variant",
                          "splice_donor_region_variant", "splice_polypyrimidine_tract_variant")

KF_TYPES = c("small variant", "SV", "CNA")

KF_SV_CLASS_LABEL = c(translocation = "translocation",
                      `intra-chr breakend` = "intrachromosomal rearrangement",
                      `single breakend` = "single breakend",
                      DEL = "deletion", DUP = "duplication", INV = "inversion",
                      INS = "insertion", BND = "breakend", sBND = "single breakend")

KF_CALLER_LABEL = c(severus = "Severus", savana = "SAVANA")

empty_findings = function() {
  data.table(gene = character(), type = character(), alteration = character(),
             evidence = character(), tier = integer(), role = character(),
             panels = character(), vaf = numeric(), source_id = character(),
             genes = character(),
             dp = numeric(), callers = character(), locus = character(),
             loh = logical(), in_panel = logical(), rescued = logical(), ev_tier = integer())
}

# ---- tiering -------------------------------------------------------------------------------

# Tier of a finding from its type, its evidence tokens and the gene's role. Pure, vectorised.
#   type      "small variant" | "SV" | "CNA"
#   evidence  "+"-joined tokens:
#               small variant: clinvar (ClinVar P/LP), lof (HIGH truncating), high (other HIGH),
#                              missense (MODERATE + in-silico damaging), biallelic (LOH in gene)
#               SV:  direct (breakpoint in the gene, or a focal DEL/DUP covering it),
#                    padfoot (Padfoot HIGH*), padfoot_low (Padfoot annotated the call and
#                    rated none of it HIGH: direct drops to tier 2), proximity (only within
#                    the window),
#                    recombining (a non-translocation event <= KF_FOCAL_SV_MAX at an IG/TR
#                    locus: V(D)J recombination, tier 3), intronic (a
#                    DEL/DUP/INS VEP annotated without exon effect, no Padfoot HIGH: tier 3),
#                    unannotated (such a call with neither VEP nor Padfoot: tier 2),
#                    rescued (every call of the event whitelist-rescued: KF_TIER_RESCUED)
#               CNA: amplification, homdel, par (pseudoautosomal), recombining (IG/TR locus)
#   role      "oncogene" | "tsg" | "both" | "" / NA
# Returns integer 1-3, NA = not listed (amplification of a tumour suppressor).
kf_tier = function(type, evidence, role) {
  n = max(length(type), length(evidence), length(role))
  if (n == 0) return(integer(0))
  type = rep_len(as.character(type), n)
  evidence = rep_len(as.character(evidence), n)
  role = rep_len(as.character(role), n)
  role[is.na(role)] = ""
  vapply(seq_len(n), function(i) {
    tok = strsplit(evidence[i], "+", fixed = TRUE)[[1]]
    has = function(t) t %in% tok
    r = role[i]
    switch(type[i],
      "small variant" = {
        t = if (has("clinvar") || has("lof") || has("high")) 1L
            else if (has("missense")) 2L else NA_integer_
        if (is.na(t)) return(NA_integer_)
        # A truncating hit in an oncogene is not how an oncogene is activated
        if (r == "oncogene" && has("lof") && !has("clinvar")) t = 2L
        # Second hit in a tumour suppressor: the missense tier is lifted
        if (r == "tsg" && has("biallelic")) t = 1L
        t
      },
      "SV" = {
        t = if (has("recombining") || has("intronic")) 3L
            else if (has("padfoot")) 1L
            # Padfoot looked at the breakpoint and did not rate it HIGH
            else if (has("direct") && (has("padfoot_low") || has("unannotated"))) 2L
            else if (has("direct")) 1L
            else if (has("proximity")) 3L else NA_integer_
        if (!is.na(t) && has("rescued")) KF_TIER_RESCUED else t
      },
      "CNA" = {
        if (has("par") || has("recombining")) return(3L)
        if (has("amplification")) {
          if (r == "tsg") NA_integer_ else 1L
        } else if (has("homdel")) {
          if (r == "oncogene") 2L else 1L
        } else NA_integer_
      },
      NA_integer_)
  }, integer(1))
}

# The note a role adjustment adds to the evidence text ("" when none applied)
kf_tier_note = function(type, evidence, role) {
  n = max(length(type), length(evidence), length(role))
  if (n == 0) return(character(0))
  type = rep_len(as.character(type), n)
  evidence = rep_len(as.character(evidence), n)
  role = rep_len(as.character(role), n)
  role[is.na(role)] = ""
  vapply(seq_len(n), function(i) {
    tok = strsplit(evidence[i], "+", fixed = TRUE)[[1]]
    if (type[i] == "small variant" && role[i] == "oncogene" && "lof" %in% tok &&
        !"clinvar" %in% tok) return("loss-of-function in an oncogene")
    if (type[i] == "CNA" && role[i] == "oncogene" && "homdel" %in% tok &&
        !any(c("par", "recombining") %in% tok)) return("deletion of an oncogene")
    ""
  }, character(1))
}

# ---- panel helpers --------------------------------------------------------------------------

# Upper-cased SNV-scope / SV-scope symbols of one panel (after the YAML round trip)
kf_panel_snv = function(p) panel_snv_genes(list(p))
kf_panel_sv = function(p) {
  g = unlist(if (!is.null(p$genes_sv)) p$genes_sv
             else if (isTRUE(as.logical(unlist(p$has_scopes)))) NULL else p$genes)
  g = toupper(trimws(as.character(g)))
  unique(g[!is.na(g) & nzchar(g)])
}

# One role per gene over every panel: data.table(gene (upper), role). Panels that disagree
# (oncogene in one, tsg in another) give "both", which applies no adjustment.
kf_gene_roles = function(all_panels) {
  r = rbindlist(c(list(data.table(gene = character(), role = character())),
                  lapply(all_panels, function(p) if (is.list(p)) panel_roles(p))))
  if (nrow(r) == 0) return(data.table(gene = character(), role = character()))
  r = r[nzchar(role) & !is.na(role)]
  r[, gene := toupper(gene)]
  r[, .(role = if (uniqueN(role) == 1) role[1] else "both"), by = gene]
}

kf_role_of = function(genes, roles) {
  out = roles$role[match(toupper(genes), roles$gene)]
  out[is.na(out)] = ""
  out
}

# Comma-joined names of the panels whose set (named list of upper-case vectors) has each gene
kf_panels_of = function(genes, sets) {
  g = toupper(genes)
  if (length(sets) == 0) return(rep("", length(g)))
  m = vapply(sets, function(s) g %in% s, logical(length(g)))
  if (!is.matrix(m)) m = matrix(m, nrow = length(g))
  apply(m, 1, function(r) paste(names(sets)[r], collapse = ","))
}

# Union of comma-joined panel strings
kf_union_panels = function(x, order = NULL) {
  p = unique(unlist(strsplit(x[!is.na(x)], ",", fixed = TRUE)))
  p = p[nzchar(p)]
  if (length(order)) p = c(intersect(order, p), setdiff(p, order))
  paste(p, collapse = ",")
}

# ---- formatting ----------------------------------------------------------------------------

# "ENSP00000301067.7:p.Arg2645Ter" -> "p.Arg2645Ter"; "%3D" (VEP's escaped "=") -> "="
kf_short_hgvsp = function(x) {
  x = as.character(x)
  x = sub("^[^:]*:(?=p\\.)", "", x, perl = TRUE)
  x = gsub("%3D", "=", x, fixed = TRUE)
  x[!is.na(x) & !nzchar(trimws(x))] = NA_character_
  x
}

kf_pretty_csq = function(x) {
  x = gsub("_variant", "", as.character(x), fixed = TRUE)
  x = gsub("[,&]+", ", ", x)
  gsub("_", " ", x, fixed = TRUE)
}

kf_int = function(x) ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "d", big.mark = ","))

# Join non-empty pieces row-wise
kf_join = function(..., sep = " · ") {
  parts = list(...)
  n = max(vapply(parts, length, integer(1)))
  parts = lapply(parts, function(p) { p = rep_len(as.character(p), n); p[is.na(p)] = ""; p })
  vapply(seq_len(n), function(i) {
    v = vapply(parts, `[`, character(1), i)
    paste(v[nzchar(v)], collapse = sep)
  }, character(1))
}

# Does [start, end] on chrom overlap a pseudoautosomal region of `reference`?
kf_in_par = function(chrom, start, end, reference) {
  ref = normalise_reference_name(reference)
  par = if (!is.na(ref)) KF_PAR[[ref]] else NULL
  n = length(chrom)
  if (is.null(par) || n == 0) return(rep(FALSE, n))
  chrom = ensure_chr_prefix(as.character(chrom))
  vapply(seq_len(n), function(i) {
    if (is.na(start[i]) || is.na(end[i])) return(FALSE)
    any(par$chrom == chrom[i] & par$start <= end[i] & par$end >= start[i])
  }, logical(1))
}

# ---- small variants ------------------------------------------------------------------------

# Evidence flags per variant_table row. Every optional annotation column is guarded: absent
# means the rule cannot fire (the plugin did not run), never an error.
kf_small_variant_flags = function(vt) {
  n = nrow(vt)
  col = function(nm) if (nm %in% names(vt)) as.character(vt[[nm]]) else rep(NA_character_, n)
  lc  = function(nm) tolower(trimws(col(nm)))
  impact = toupper(col("impact"))
  csq    = col("consequence")
  clin   = if ("clinvar" %in% names(vt)) clinvar_is_pathogenic(vt$clinvar) else rep(FALSE, n)
  high   = !is.na(impact) & impact == "HIGH"
  mod    = !is.na(impact) & impact == "MODERATE"
  terms  = strsplit(ifelse(is.na(csq), "", csq), "[,&]+")
  # ClinVar classifies the variant, not each gene VEP lists for it: when some of a variant's
  # rows are coding/splice, only those carry it (ATM p.Arg23Ter, not the Y_RNA downstream)
  if (any(clin) && all(c("chrom", "pos", "ref", "alt") %in% names(vt))) {
    in_gene = vapply(terms, function(t) any(trimws(t) %in% KF_CLINVAR_GENE_TERMS), logical(1))
    key = paste(vt$chrom, vt$pos, vt$ref, vt$alt)
    any_in_gene = ave(in_gene, key, FUN = any)
    clin = clin & (in_gene | !any_in_gene)
  }
  lof    = high & vapply(terms, function(t) any(trimws(t) %in% KF_LOF_TERMS), logical(1))
  am     = lc("am_class")
  am_lp  = !is.na(am) & am == "likely_pathogenic"
  # AlphaMissense first; a row it does not cover (column absent, or no call for this row)
  # falls back to EVE, then to SIFT and PolyPhen agreeing
  eve    = is.na(am) & !is.na(lc("eve_class")) & lc("eve_class") %in% c("pathogenic", "likely_pathogenic")
  sp     = is.na(am) & !eve & !is.na(lc("sift")) & lc("sift") == "deleterious" &
           !is.na(lc("polyphen")) & lc("polyphen") == "probably_damaging"
  missense = mod & (am_lp | eve | sp)
  list(clin = clin, high = high, lof = lof, missense = missense,
       am_lp = mod & am_lp, eve = mod & eve, sp = mod & sp)
}

kf_small_variant_rows = function(vt, snv_sets, roles, loh_genes, max_nonpanel = KF_MAX_NONPANEL) {
  out = empty_findings()
  attr(out, "n_nonpanel") = list(kept = 0L, total = 0L)
  if (is.null(vt) || nrow(vt) == 0 || !"symbol" %in% names(vt)) return(out)
  vt = as.data.table(vt)
  sym = toupper(trimws(as.character(vt$symbol)))
  ok_sym = !is.na(sym) & nzchar(sym) & sym != "-"
  f = kf_small_variant_flags(vt)
  panel_genes = unique(unlist(snv_sets, use.names = FALSE))
  in_panel = ok_sym & sym %in% panel_genes
  keep_panel = in_panel & (f$clin | f$high | f$missense)
  keep_other = ok_sym & !in_panel & (f$clin | f$high)

  total_other = sum(keep_other)
  idx_other = which(keep_other)
  if (length(idx_other) > max_nonpanel) {
    vaf = if ("vaf" %in% names(vt)) suppressWarnings(as.numeric(vt$vaf)) else rep(NA_real_, nrow(vt))
    o = order(-(f$clin[idx_other] & f$high[idx_other]), -f$clin[idx_other],
              -f$high[idx_other], -ifelse(is.na(vaf[idx_other]), -1, vaf[idx_other]))
    idx_other = idx_other[o][seq_len(max_nonpanel)]
  }
  idx = c(which(keep_panel), idx_other)
  if (length(idx) == 0) {
    attr(out, "n_nonpanel") = list(kept = 0L, total = as.integer(total_other))
    return(out)
  }

  d = vt[idx]
  g = sym[idx]
  ff = lapply(f, `[`, idx)
  col = function(nm, na = NA_character_) if (nm %in% names(d)) d[[nm]] else rep(na, nrow(d))
  csq = as.character(col("consequence"))
  hg = kf_short_hgvsp(col("hgvsp"))
  ref = as.character(col("ref")); alt = as.character(col("alt"))
  locus = paste0(col("chrom"), ":", kf_int(col("pos")))
  alteration = ifelse(!is.na(hg), hg,
                      paste0(kf_pretty_csq(csq), " ", ref, ">", alt))

  am_score = suppressWarnings(as.numeric(col("am_score", NA_real_)))
  loh = g %in% loh_genes
  biallelic = loh & (in_panel[idx])
  ev_clin = ifelse(ff$clin, paste0("ClinVar: ", gsub("_", " ", col("clinvar"))), "")
  ev_high = ifelse(ff$high, paste0("HIGH impact (", kf_pretty_csq(csq), ")"), "")
  ev_mis  = ifelse(ff$am_lp,
                   paste0("AlphaMissense likely pathogenic",
                          ifelse(is.na(am_score), "", sprintf(" (%.2f)", am_score))),
            ifelse(ff$eve, "EVE pathogenic",
            ifelse(ff$sp, "SIFT deleterious + PolyPhen probably damaging", "")))
  ev_mis  = ifelse(nzchar(ev_mis), paste0("MODERATE impact, ", ev_mis), "")
  ev_bi   = ifelse(biallelic, "possible biallelic (LOH)", "")

  tokens = kf_join(ifelse(ff$clin, "clinvar", ""),
                   ifelse(ff$lof, "lof", ifelse(ff$high, "high", "")),
                   ifelse(ff$missense & !ff$high, "missense", ""),
                   ifelse(biallelic, "biallelic", ""), sep = "+")
  role = ifelse(in_panel[idx], kf_role_of(g, roles), "")
  tier = kf_tier("small variant", tokens, role)
  note = kf_tier_note("small variant", tokens, role)

  out = data.table(
    gene = g, type = "small variant", alteration = alteration,
    evidence = kf_join(ev_clin, ev_high, ev_mis, ev_bi, note),
    tier = as.integer(tier), role = role,
    panels = kf_panels_of(g, snv_sets),
    vaf = suppressWarnings(as.numeric(col("vaf", NA_real_))),
    source_id = paste0(col("chrom"), ":", col("pos"), ":", ref, ">", alt),
    genes = g,
    dp =suppressWarnings(as.numeric(col("dp", NA_real_))),
    callers = as.character(col("callers", "")),
    locus = locus, loh = loh, in_panel = in_panel[idx])
  out = out[!is.na(tier)]
  attr(out, "n_nonpanel") = list(kept = length(idx_other), total = as.integer(total_other))
  out
}

# ---- structural variants -------------------------------------------------------------------

# Symbols of a Padfoot gene cell: "EWSR1 (intron7), ATF1 (intron3)" -> c("EWSR1", "ATF1")
kf_padfoot_symbols = function(x) {
  lapply(as.character(x), function(cell) {
    if (is.na(cell) || !nzchar(cell)) return(character(0))
    s = trimws(sub("\\s*\\(.*$", "", trimws(strsplit(cell, ",(?![^(]*\\))", perl = TRUE)[[1]])))
    unique(toupper(s[nzchar(s)]))
  })
}

# One short locus string per SV row: "chr12–chr22" for a translocation, else chrom:start–end Mb
kf_sv_locus = function(t) {
  ca = as.character(t$chrom_a); cb = as.character(t$chrom_b)
  pa = suppressWarnings(as.numeric(t$pos_a)); pb = suppressWarnings(as.numeric(t$pos_b))
  cls = as.character(t$svclass)
  rk = function(ch) { r = match(ch, chromosomes_for_sex("male")); r[is.na(r)] = 999L; r }
  ifelse(cls == "translocation",
         ifelse(rk(ca) <= rk(cb), paste0(ca, "–", cb), paste0(cb, "–", ca)),
  ifelse(cls == "single breakend" | is.na(pb),
         paste0(ca, ":", sprintf("%.2f", pa / 1e6), " Mb"),
         paste0(ca, ":", sprintf("%.2f", pmin(pa, pb) / 1e6), "–",
                sprintf("%.2f", pmax(pa, pb) / 1e6), " Mb")))
}

# Is either breakpoint of each hit inside an interval of its gene? `iv` = data.table(chrom,
# start, end, gene (upper)) over every panel
kf_bp_in_gene = function(chrom_a, pos_a, chrom_b, pos_b, gene, iv) {
  n = length(gene)
  if (n == 0 || is.null(iv) || nrow(iv) == 0) return(rep(FALSE, n))
  vapply(seq_len(n), function(i) {
    want = gene[i]   # outside `[`: inside it `gene` and `i` would mean iv's own column / row
    g = iv[iv$gene == want]
    if (nrow(g) == 0) return(FALSE)
    inside = function(ch, p) !is.na(ch) && !is.na(p) && any(g$chrom == ch & g$start <= p & g$end >= p)
    inside(chrom_a[i], pos_a[i]) || inside(chrom_b[i], pos_b[i])
  }, logical(1))
}

# Does the span [min(pos_a, pos_b), max(pos_a, pos_b)] of a same-contig SV overlap `gene`'s
# interval in `iv`? Vectorised over rows like kf_bp_in_gene()
kf_span_covers_gene = function(chrom_a, pos_a, chrom_b, pos_b, gene, iv) {
  n = length(gene)
  if (n == 0 || is.null(iv) || nrow(iv) == 0) return(rep(FALSE, n))
  vapply(seq_len(n), function(i) {
    want = gene[i]
    ch = chrom_a[i]
    if (is.na(ch) || is.na(chrom_b[i]) || ch != chrom_b[i] || is.na(pos_a[i]) || is.na(pos_b[i]))
      return(FALSE)
    g = iv[iv$gene == want]
    any(g$chrom == ch & g$start <= max(pos_a[i], pos_b[i]) & g$end >= min(pos_a[i], pos_b[i]))
  }, logical(1))
}

# SV classes that say only "a junction here", not what kind of event it is
KF_SV_UNTYPED = c("intra-chr breakend", "single breakend", "BND", "sBND")

# Which rows of an SV table are the same event called by different callers.
# Severus and SAVANA report the same rearrangement with breakpoints a few bp apart, and type it
# differently: SAVANA's same-chromosome pairs are typed DEL/DUP/INV from BP_NOTATION, Severus
# leaves many as an untyped "intra-chr breakend". Rows are matched when they come from
# different callers, join the same contig(s), and both breakpoints agree within `tol` (on one
# contig also within half the smaller event, so two neighbouring small deletions stay apart).
# Matches are transitive (union-find). Returns data.table(event (integer id, the smallest row
# index of the event), ev_class): an untyped row takes its event's typed class when the event
# has exactly one.
kf_sv_events = function(st, tol = KF_SV_SAME_EVENT_BP) {
  n = nrow(st)
  if (n == 0) return(data.table(event = integer(0), ev_class = character(0)))
  chroms = chromosomes_for_sex("male")
  rk = function(ch) { r = match(ch, chroms); r[is.na(r)] = 999L; r }
  ca = as.character(st$chrom_a); cb = as.character(st$chrom_b)
  pa = suppressWarnings(as.numeric(st$pos_a)); pb = suppressWarnings(as.numeric(st$pos_b))
  cls = as.character(st$svclass)
  caller = if ("caller" %in% names(st)) as.character(st$caller) else rep("", n)
  caller[is.na(caller)] = ""
  single = cls %in% c("single breakend", "sBND") | is.na(cb) | is.na(pb)
  # Orient every pair the same way: lower-ranked contig first, then lower position
  swap = !single & !is.na(ca) &
         (rk(cb) < rk(ca) | (rk(cb) == rk(ca) & cb < ca) | (cb == ca & pb < pa))
  swap[is.na(swap)] = FALSE
  c1 = ifelse(swap, cb, ca); p1 = ifelse(swap, pb, pa)
  c2 = ifelse(swap, ca, cb); p2 = ifelse(swap, pa, pb)
  c2[single] = NA_character_; p2[single] = NA_real_

  parent = seq_len(n)
  root = function(i) { while (parent[i] != i) i = parent[i]; i }
  for (k in split(seq_len(n), paste(c1, c2, single))) {
    if (length(k) < 2) next
    for (a in seq_len(length(k) - 1L)) for (b in (a + 1L):length(k)) {
      i = k[a]; j = k[b]
      if (caller[i] == caller[j] || is.na(p1[i]) || is.na(p1[j])) next
      tl = tol
      if (!single[i] && identical(c1[i], c2[i]))
        tl = min(tl, max(10, 0.5 * min(abs(p2[i] - p1[i]), abs(p2[j] - p1[j]))))
      if (abs(p1[i] - p1[j]) > tl) next
      if (!single[i] && (is.na(p2[i]) || is.na(p2[j]) || abs(p2[i] - p2[j]) > tl)) next
      ri = root(i); rj = root(j)
      if (ri != rj) parent[max(ri, rj)] = min(ri, rj)
    }
  }
  event = vapply(seq_len(n), root, integer(1))

  ev_class = cls
  for (e in unique(event[duplicated(event)])) {
    m = which(event == e)
    typed = unique(cls[m][!is.na(cls[m]) & !cls[m] %in% KF_SV_UNTYPED])
    if (length(typed) == 1) ev_class[m] = typed
  }
  data.table(event = event, ev_class = ev_class)
}

kf_sv_rows = function(sv_table, all_panels, roles, max_proximity = KF_MAX_PROXIMITY_PER_GENE) {
  out = empty_findings()
  if (is.null(sv_table) || nrow(sv_table) == 0 || length(all_panels) == 0) return(out)
  st = copy(as.data.table(sv_table))
  need = c("chrom_a", "pos_a", "chrom_b", "pos_b", "svtype")
  if (!all(need %in% names(st))) return(out)
  if (!"svclass" %in% names(st)) st[, svclass := svtype]
  if (!"caller" %in% names(st)) st[, caller := "severus"]
  if (!"id" %in% names(st)) st[, id := paste0("sv", .I)]
  st[, `:=`(pos_a = suppressWarnings(as.numeric(pos_a)), pos_b = suppressWarnings(as.numeric(pos_b)))]
  ev = kf_sv_events(st)
  st[, `:=`(kf_event = ev$event, kf_class = ev$ev_class)]
  # An event is rescued only when every call of it is: a rescued Severus call that SAVANA also
  # made stands with SAVANA's
  st[, kf_rescued := sv_is_rescued(st)]
  st[, kf_ev_rescued := all(kf_rescued), by = kf_event]

  # Coordinate / symbol hits, every panel at once (SV scope is applied inside)
  ht = sv_panel_hit_table(st, all_panels)
  hits = if (nrow(ht) > 0)
    ht[, .(row, gene = toupper(gene), panel, direct, gap, side, match, padfoot = FALSE)]
  else data.table(row = integer(), gene = character(), panel = character(),
                  direct = logical(), gap = numeric(), side = character(),
                  match = character(), padfoot = logical())

  # sv_panel_hit_table() calls a DEL/DUP/INV/INS "direct" whenever its span overlaps the gene;
  # a 200 Mb inversion then "hits" every panel gene on the arm. Here a span hit counts only
  # when a breakpoint is in the gene, or the event is a focal (<= KF_FOCAL_SV_MAX) DEL/DUP
  # covering it; a gene merely inside a large event is left to the copy-number rows.
  # A symbol hit is the same case: VEP lists every gene a record spans (a Severus
  # same-chromosome BND included), so it is position-tested too wherever a list gives the
  # gene coordinates
  hits[, focal := FALSE]
  iv = rbindlist(lapply(all_panels, function(p) panel_gene_intervals(p, "sv")))
  if (nrow(iv)) iv[, gene := toupper(gene)]
  sp = which(hits$direct %in% TRUE &
             (hits$side == "span" | (hits$match == "symbol" & hits$gene %in% iv$gene)))
  if (length(sp)) {
    r = hits$row[sp]
    # (no local named like a column: data.table's i/j would read the column instead)
    is_bp = kf_bp_in_gene(st$chrom_a[r], st$pos_a[r], st$chrom_b[r], st$pos_b[r], hits$gene[sp], iv)
    sz = abs(st$pos_b[r] - st$pos_a[r])
    # A span hit overlaps the gene by construction; a VEP symbol may name a gene up to 5 kb away
    covers = hits$side[sp] == "span" |
             kf_span_covers_gene(st$chrom_a[r], st$pos_a[r], st$chrom_b[r], st$pos_b[r], hits$gene[sp], iv)
    # The event's class: a Severus breakend matched to SAVANA's DEL is a deletion
    is_focal = !is_bp & covers & st$kf_class[r] %in% c("DEL", "DUP") &
               !is.na(sz) & sz <= KF_FOCAL_SV_MAX
    set(hits, i = sp, j = "focal", value = is_focal)
    drop_i = sp[!is_bp & !is_focal]
    if (length(drop_i)) set(hits, i = drop_i, j = "direct", value = NA)   # dropped below
  }
  hits = hits[!is.na(direct)]

  # Padfoot: a HIGH* call whose annotated genes include a panel gene counts as strong evidence
  if ("padfoot_impact" %in% names(st) && "padfoot_genes" %in% names(st)) {
    hi = which(!is.na(st$padfoot_impact) & startsWith(toupper(st$padfoot_impact), "HIGH"))
    if (length(hi)) {
      pg = kf_padfoot_symbols(st$padfoot_genes[hi])
      sets = lapply(all_panels, kf_panel_sv)
      pf = rbindlist(lapply(seq_along(hi), function(k) {
        rbindlist(lapply(names(sets), function(nm) {
          gg = intersect(pg[[k]], sets[[nm]])
          if (!length(gg)) return(NULL)
          data.table(row = hi[k], gene = gg, panel = nm, direct = FALSE, gap = NA_real_,
                     side = "", match = "padfoot", padfoot = TRUE, focal = FALSE)
        }))
      }))
      if (nrow(pf)) hits = rbind(hits, pf)
    }
  }
  if (nrow(hits) == 0) return(out)

  # One record per (SV row, gene): strongest evidence, nearest gap, every panel
  rg = hits[, .(bp = any(direct & !focal, na.rm = TRUE), focal = any(direct & focal, na.rm = TRUE),
                padfoot = any(padfoot),
                gap = suppressWarnings(min(gap, na.rm = TRUE)),
                side = { s = side[nzchar(side)]; if (length(s)) s[1] else "A" },
                panels = paste(unique(panel), collapse = ",")), by = .(row, gene)]
  rg[!is.finite(gap), gap := NA_real_]
  rg[, strong := bp | focal | padfoot]
  r = rg$row
  col = function(nm, na = NA_character_)
    if (nm %in% names(st)) as.character(st[[nm]])[r] else rep(na, length(r))
  # The event's class: a Severus "intra-chr breakend" matched to SAVANA's DEL reads as a deletion
  rg[, svclass := col("kf_class")]
  rg[, event := st$kf_event[r]]
  # Callers label one junction differently (Severus "intra-chr breakend", SAVANA "INV"): one group
  rg[, grp := fifelse(svclass %in% c("intra-chr breakend", "INV"), "intrachromosomal", svclass)]
  rg[, caller := col("caller")]
  rg[, id := col("id")]
  rg[, vaf := suppressWarnings(as.numeric(col("vaf")))]
  rg[, locus := kf_sv_locus(st[r])]
  rg[, size := abs(st$pos_b[r] - st$pos_a[r])]
  rg[, gchrom := fifelse(side == "B", as.character(st$chrom_b[r]), as.character(st$chrom_a[r]))]
  rg[, pf_impact := col("padfoot_impact")]
  rg[, pf_genes := col("padfoot_genes")]
  rg[, vep_impact := toupper(col("impact"))]
  rg[, vep_csq := col("consequence")]
  rg[, rescued := st$kf_ev_rescued[r]]
  rg[, support := suppressWarnings(as.integer(col("support")))]
  rg[, filt := col("filter")]
  rg[, wl := col("wl_rescue")]
  # Partner genes of a junction: Padfoot's, else VEP's symbols on either side
  pf_sym = kf_padfoot_symbols(rg$pf_genes)
  vep_sym = kf_padfoot_symbols(paste(col("gene_a", ""), col("gene_b", ""), sep = ","))
  rg[, partners := mapply(function(a, b, g) {
    s = if (length(a)) a else b
    s = s[!is.na(s) & s != "NA"]
    paste(setdiff(s, g), collapse = ", ")
  }, pf_sym, vep_sym, gene)]

  chroms = chromosomes_for_sex("male")
  class_lab = function(x) {
    x = unique(x)
    l = KF_SV_CLASS_LABEL[x]; l[is.na(l)] = x[is.na(l)]
    if (all(x %in% c("intra-chr breakend", "INV")) && length(x) > 1)
      "intrachromosomal rearrangement / inversion"
    else paste(unique(unname(l)), collapse = " / ")
  }

  # The pieces of each (gene, event group) row; the text is composed after a rearrangement
  # joining two listed genes has been merged into one row
  agg = rg[, {
    ord = order(-(bp | padfoot), gap, na.last = TRUE)
    loci = unique(locus[ord])
    junction = grp[1] %in% c("translocation", "intrachromosomal")
    prt = unique(unlist(strsplit(partners[nzchar(partners)], ", ", fixed = TRUE)))
    # Padfoot annotated these calls and rated none HIGH: an intronic / benign-looking breakpoint
    # (the call itself rated low: a call Padfoot rated HIGH for its other gene is not low here)
    pf_low = !any(padfoot) && any(!is.na(pf_impact) & !startsWith(toupper(pf_impact), "HIGH"))
    # A small DEL/DUP/INS that VEP annotated as touching no exon (impact below HIGH on every
    # call; mostly 50-300 bp intronic deletions) is not a disruption. Junctions are exempt:
    # an intronic breakpoint is exactly how a fusion looks
    vep_low = !any(padfoot) && !junction && !any(focal & !bp) &&
              all(svclass %in% c("DEL", "DUP", "INS")) &&
              any(!is.na(vep_impact)) && !any(vep_impact == "HIGH", na.rm = TRUE)
    # Neither VEP nor Padfoot annotated such a call: whether it touches an exon is unknown
    unannot = !any(padfoot) && !junction && !any(focal & !bp) &&
              all(svclass %in% c("DEL", "DUP", "INS")) &&
              all(is.na(vep_impact)) && all(is.na(pf_impact))
    main = if (!strong[1]) paste0("near ", gene[1], ", ", fmt_bp(min(gap, na.rm = TRUE)), " away")
           else if (any(bp)) paste0("breakpoint in ", gene[1])
           else if (any(focal)) paste0("focal ", class_lab(svclass[focal]), " (",
                                       fmt_bp(min(size[focal], na.rm = TRUE)), ") covers ", gene[1])
           else ""
    pf_txt = if (!strong[1]) ""
             else if (any(padfoot)) {
               k = which(padfoot)[1]
               paste0("Padfoot ", pf_impact[k],
                      if (!is.na(pf_genes[k]) && nzchar(pf_genes[k])) paste0(": ", pf_genes[k]) else "")
             } else if (pf_low) {
               paste0("Padfoot rates it ", paste(sort(unique(pf_impact[!is.na(pf_impact)])), collapse = "/"))
             } else ""
    low_txt = if (!strong[1]) ""
              else if (vep_low) paste0("no exon affected (VEP: ",
                                       kf_pretty_csq(paste(unique(vep_csq[!is.na(vep_csq)]), collapse = ",")), ")")
              else if (unannot) "exon effect not annotated (no VEP / Padfoot annotation for this call)"
              else ""
    rk = match(gchrom, chroms); rk[is.na(rk)] = 999L
    resc_txt = if (!rescued) "" else {
      w = sort(unique(wl[!is.na(wl)]))
      failed = !is.na(filt) & !filt %in% SV_PASS_FILTERS
      f = sort(unique(filt[failed]))
      s = support[!is.na(support)]
      n = uniqueN(id); nf = uniqueN(id[failed])
      paste0(if (n == 1L) "whitelist-rescued" else paste(n, "calls whitelist-rescued"),
             if (length(w)) paste0(" (", paste(w, collapse = "/"), ")"),
             if (nf > 0L) paste0(", ", if (n == 1L) "" else paste0(nf, " "), "failed ",
                                 paste(f, collapse = "/")),
             if (length(s)) paste0(" · ", if (min(s) == max(s)) min(s) else paste0(min(s), "–", max(s)),
                                   " supporting read", if (max(s) == 1L) "" else "s",
                                   if (n > 1L && min(s) == max(s)) " each" else ""))
    }
    list(cls = class_lab(svclass),
         loc = paste0(loci[1], if (length(loci) > 1) paste0(" +", length(loci) - 1L, " more") else ""),
         partner = paste(prt, collapse = ", "),
         callers = paste(unique(caller), collapse = ","),
         junction = junction, any_bp = any(bp), main = main, pf_txt = pf_txt, low_txt = low_txt,
         resc_txt = resc_txt,
         pf_first = if (any(padfoot)) pf_genes[which(padfoot)[1]] else NA_character_,
         crank = min(rk),
         events = paste(sort(unique(event)), collapse = ","),
         tokens = if (strong[1])
                    paste(c(if (any(bp | focal)) "direct", if (any(padfoot)) "padfoot",
                            if (pf_low) "padfoot_low", if (vep_low) "intronic",
                            if (unannot) "unannotated",
                            # V(D)J recombination: any non-translocation event of locus size
                            if (grepl(KF_RECOMBINING_LOCI, gene[1]) && grp[1] != "translocation" &&
                                all(!is.na(size) & size <= KF_FOCAL_SV_MAX)) "recombining",
                            if (rescued) "rescued"),
                          collapse = "+")
                  else paste(c("proximity", if (rescued) "rescued"), collapse = "+"),
         panels = kf_union_panels(panels, names(all_panels)),
         vaf = suppressWarnings(max(vaf, na.rm = TRUE)),
         source_id = paste(unique(id), collapse = ","),
         locus = loci[1],
         min_gap = suppressWarnings(min(gap, na.rm = TRUE)))
  }, by = .(gene, grp, strong, rescued)]
  agg[!is.finite(vaf), vaf := NA_real_]

  # A gene's proximity rows: when the same event (any caller) already hits it strongly they
  # add nothing, and a noisy sample is capped per gene (nearest first)
  ev_of = function(x) strsplit(x, ",", fixed = TRUE)[[1]]
  strong_ev = agg[strong == TRUE, .(ev = paste(events, collapse = ",")), by = gene]
  redundant = mapply(function(g, e) {
    s = strong_ev$ev[strong_ev$gene == g]
    length(s) > 0 && all(ev_of(e) %in% ev_of(s))
  }, agg$gene, agg$events)
  # A rescued call merely near a gene would be tier 3: never listed, so never counted in the cap
  agg = agg[(strong == TRUE | !redundant) & !(rescued & !strong)]
  agg[, n_more := 0L]
  prox = agg[strong == FALSE][order(gene, min_gap)]
  if (nrow(prox) > 0) {
    prox[, k := seq_len(.N), by = gene]
    prox[, n_g := .N, by = gene]
    prox[k == max_proximity & n_g > max_proximity, n_more := n_g - max_proximity]
    prox = prox[k <= max_proximity][, c("k", "n_g") := NULL]
  }
  agg = rbind(agg[strong == TRUE], prox)
  if (nrow(agg) == 0) return(out)

  agg[, role := kf_role_of(gene, roles)]
  agg[, tier := kf_tier("SV", tokens, role)]
  agg[, ev_tier := kf_tier("SV", gsub("+rescued", "", tokens, fixed = TRUE), role)]
  # A rescued call that would be tier 3 on its evidence is not listed at all
  agg = agg[!is.na(tier) & !(rescued & ev_tier >= 3L)]
  if (nrow(agg) == 0) return(out)
  agg[, `:=`(genes = gene, merged = FALSE)]
  agg = kf_merge_sv_partners(agg, names(all_panels))

  # Compose the text. Callers are those of the row's events, not only of the calls that hit
  # the gene: Severus' copy of a SAVANA deletion may only reach the gene as "near"
  caller_lab = function(x) { l = KF_CALLER_LABEL[x]; ifelse(is.na(l), x, l) }
  ev_callers = split(as.character(st$caller), st$kf_event)
  cl = mapply(function(x, e) {
    x = c(strsplit(x, ",", fixed = TRUE)[[1]],
          unlist(ev_callers[strsplit(e, ",", fixed = TRUE)[[1]]], use.names = FALSE))
    x = sort(unique(x[!is.na(x) & nzchar(x)]))
    paste(unique(caller_lab(x)), collapse = ", ")
  }, agg$callers, agg$events, USE.NAMES = FALSE)
  n_ev = lengths(strsplit(agg$events, ",", fixed = TRUE))
  show_partner = agg$strong & agg$junction & !agg$merged & nzchar(agg$partner)
  alteration = paste0(agg$cls, " ", agg$loc,
                      ifelse(show_partner,
                             paste0(", partner ", vapply(strsplit(agg$partner, ", ", fixed = TRUE),
                                                         function(p) paste(utils::head(p, 3), collapse = "/"),
                                                         character(1))), ""),
                      " (", cl, ")")
  evidence = kf_join(agg$resc_txt, agg$main, agg$pf_txt, agg$low_txt)
  evidence = paste0(evidence, ifelse(n_ev > 1, paste0(" · ", n_ev, " events"), ""))
  evidence = paste0(evidence, ifelse(agg$n_more > 0,
                                     paste0(" (+", agg$n_more, " more nearby SV group",
                                            ifelse(agg$n_more == 1, "", "s"), " not listed)"), ""))
  # A locus-sized non-translocation event at IG/TR is V(D)J recombination, not a driver
  evidence = ifelse(grepl("recombining", agg$tokens),
                    paste0(evidence, " · within an immunoglobulin / T-cell receptor locus — physiological rearrangement"),
                    evidence)

  data.table(gene = agg$gene, type = "SV", alteration = alteration,
             evidence = evidence, tier = as.integer(agg$tier), role = agg$role,
             panels = agg$panels, vaf = agg$vaf, source_id = agg$source_id,
             genes = agg$genes, dp = NA_real_, callers = cl, locus = agg$locus,
             loh = FALSE, in_panel = TRUE, rescued = agg$rescued,
             ev_tier = as.integer(agg$ev_tier))
}

# A rearrangement whose breakpoints hit two listed genes (EWSR1-ATF1) is one row, not one per
# gene. Strong junction rows (translocation / intrachromosomal) of different genes in the same
# event group that share an event are merged, transitively. The merged row's `gene` is the
# members joined by an en dash (orientation is not established, so not HGNC "::"), in
# Padfoot's listing order where Padfoot names them, else by chromosome; `genes` is them
# comma-joined. Tier is the best member's, panels and events the union.
kf_merge_sv_partners = function(agg, panel_order = NULL) {
  j = which(agg$strong & agg$junction)
  if (length(j) < 2) return(agg)
  evs = lapply(agg$events[j], function(x) strsplit(x, ",", fixed = TRUE)[[1]])
  comp = seq_along(j)
  for (a in seq_along(j)[-length(j)]) for (b in (a + 1L):length(j)) {
    if (agg$grp[j[a]] != agg$grp[j[b]] || agg$gene[j[a]] == agg$gene[j[b]]) next
    if (length(intersect(evs[[a]], evs[[b]])) == 0) next
    old = comp[b]; comp[comp == old] = comp[a]
  }
  drop = integer(0); add = list()
  for (cc in unique(comp)) {
    m = j[comp == cc]
    if (length(unique(agg$gene[m])) < 2) next
    d = agg[m]
    pf = kf_padfoot_symbols(d$pf_first[!is.na(d$pf_first)][1])[[1]]
    o = order(match(d$gene, pf), d$crank, d$gene, na.last = TRUE)
    d = d[o]
    g = unique(d$gene)
    bp = unique(d$gene[d$any_bp])
    split_u = function(x) unique(unlist(strsplit(x[nzchar(x)], ",", fixed = TRUE)))
    tok = unique(unlist(strsplit(d$tokens, "+", fixed = TRUE)))
    add[[length(add) + 1L]] = data.table(
      gene = paste(g, collapse = "–"), grp = d$grp[1], strong = TRUE, rescued = all(d$rescued),
      cls = d$cls[1], loc = d$loc[1], partner = "",
      callers = paste(split_u(d$callers), collapse = ","),
      junction = TRUE, any_bp = length(bp) > 0,
      main = if (length(bp) == 0) ""
             else if (length(bp) == 1) paste0("breakpoint in ", bp)
             else paste0("breakpoints in ", paste(utils::head(bp, -1), collapse = ", "),
                         " and ", utils::tail(bp, 1)),
      pf_txt = { p = d$pf_txt[nzchar(d$pf_txt)]; if (length(p)) p[1] else "" },
      low_txt = "", resc_txt = { p = d$resc_txt[nzchar(d$resc_txt)]; if (length(p)) p[1] else "" },
      pf_first = d$pf_first[1], crank = min(d$crank),
      events = paste(sort(as.integer(split_u(d$events))), collapse = ","),
      tokens = paste(tok, collapse = "+"),
      panels = kf_union_panels(d$panels, panel_order),
      vaf = if (all(is.na(d$vaf))) NA_real_ else max(d$vaf, na.rm = TRUE),
      source_id = paste(split_u(d$source_id), collapse = ","),
      locus = d$locus[1], min_gap = NA_real_, n_more = 0L,
      role = "", tier = min(d$tier), ev_tier = min(d$ev_tier),
      genes = paste(g, collapse = ","), merged = TRUE)
    drop = c(drop, m)
  }
  if (length(add) == 0) return(agg)
  rbind(agg[-drop], rbindlist(add), use.names = TRUE)
}

# ---- copy number ---------------------------------------------------------------------------

kf_cna_rows = function(cna, roles, reference) {
  out = empty_findings()
  g = if (!is.null(cna)) cna$genes else NULL
  if (is.null(g) || nrow(g) == 0 || !"call" %in% names(g)) return(out)
  g = as.data.table(g)
  g = g[!is.na(call) & call %in% c("amplification", "homozygous deletion", "deep deletion")]
  if (nrow(g) == 0) return(out)
  src_lab = c(wakhan = "Wakhan", ascat = "ASCAT", savana = "SAVANA")
  if (!"span_from" %in% names(g)) g[, span_from := NA_character_]
  # One row per gene locus; a gene in several panels lists them all
  # (panels may carry slightly different spans for one gene: one row per gene and call)
  g = g[, .(panels = paste(unique(panel), collapse = ","),
            chrom = chrom[1], start = min(start), end = max(end),
            partial = any(partial %in% TRUE), loh = any(loh %in% TRUE),
            min_cn = suppressWarnings(min(min_cn, na.rm = TRUE)),
            max_cn = suppressWarnings(max(max_cn, na.rm = TRUE)),
            source = source[1], segment = any(span_from %in% "Wakhan segment")),
        by = .(gene, call)]
  par = kf_in_par(g$chrom, g$start, g$end, reference)
  recomb = grepl(KF_RECOMBINING_LOCI, toupper(g$gene))
  is_amp = g$call == "amplification"
  tokens = kf_join(ifelse(is_amp, "amplification", "homdel"),
                   ifelse(par, "par", ""), ifelse(recomb & !par, "recombining", ""), sep = "+")
  role = kf_role_of(g$gene, roles)
  tier = kf_tier("CNA", tokens, role)
  note = kf_tier_note("CNA", tokens, role)
  sl = src_lab[as.character(g$source)]
  sl[is.na(sl)] = as.character(g$source)[is.na(sl)]
  cn = ifelse(is_amp, g$max_cn, g$min_cn)
  alteration = paste0(g$call, " (CN ", trimws(formatC(cn, format = "fg", digits = 3)),
                      ifelse(g$partial, ", partial gene", ""), ")")
  evidence = kf_join(
    paste0(sl, " copy number"),
    ifelse(g$segment, "locus is Wakhan's segment, not the gene's span", ""),
    ifelse(par, "pseudoautosomal region — copy number unreliable", ""),
    ifelse(recomb & !par, "immunoglobulin / T-cell receptor locus — physiological rearrangement", ""),
    ifelse(g$loh & is_amp, "LOH", ""),
    note)
  data.table(gene = toupper(g$gene), type = "CNA", alteration = alteration,
             evidence = evidence, tier = as.integer(tier), role = role,
             panels = g$panels, vaf = NA_real_,
             source_id = paste0(g$chrom, ":", g$start, "-", g$end),
             genes = toupper(g$gene), dp = NA_real_, callers = ifelse(is.na(sl), "", sl),
             locus = paste0(ifelse(g$segment, "segment ", ""), g$chrom, ":", kf_int(g$start), "–", kf_int(g$end)),
             loh = g$loh, in_panel = TRUE)[!is.na(tier)]
}

# ---- HMF driver annotation -----------------------------------------------------------------

# The gene list built from --hmf-driver-panel carries `hmf = TRUE` (bin/render_report.R). Its
# key findings are HMF's own calls (build_hmf(), R/hmf.R) rather than the generic rules above,
# which would list each of its genes a second time under different criteria:
#   - small variants: every reported one (a reportable coding effect on the HMF canonical
#     transcript, or a hotspot), tier 1
#   - copy number: AMP and DEL from a copy-number caller tier 1, PARTIAL_AMP tier 2; a DEL called
#     only by an SV caller tier 2 up to KF_FOCAL_SV_MAX, tier 3 above it (one 70 Mb deletion
#     removes dozens of driver genes and may be heterozygous)
#   - fusions: every reported one, tier 1 at HMF likelihood HIGH (known pairs, exon
#     deletions/duplications, enhancer known pairs, known exons), tier 2 otherwise
# One row per gene and call (copy number) or per fusion pair and type, merging sources and
# callers. A fusion of two genes is "FIVE–THREE" (5′ first), the same en-dash form
# kf_merge_sv_partners() gives a rearrangement joining two listed genes.
kf_is_hmf_panel = function(p) isTRUE(as.logical(unlist(p$hmf)))

kf_hmf_rows = function(hmf, panel_names) {
  out = empty_findings()
  if (is.null(hmf) || !is.null(hmf$error) || length(panel_names) == 0) return(out)
  tag = paste(panel_names, collapse = ",")
  rows = list(out)
  na0 = function(x) { x = as.character(x); x[is.na(x)] = ""; x }

  sm = hmf$small
  if (!is.null(sm) && nrow(sm) > 0) {
    s = sm[reported %in% TRUE]
    if (nrow(s) > 0) {
      hg = kf_short_hgvsp(s$hgvsp)
      rows$small = data.table(
        gene = toupper(s$gene), type = "small variant",
        alteration = ifelse(is.na(hg), paste0(kf_pretty_csq(s$consequence), " ", s$ref, ">", s$alt), hg),
        evidence = kf_join(paste0("HMF: ", s$reason),
                           ifelse(nzchar(na0(s$hotspot)), paste0("hotspot (", s$hotspot, " match)"), ""),
                           ifelse(s$in_panel_region %in% TRUE, "ActionableCodingPanel region", "")),
        tier = 1L, role = na0(s$role), panels = tag, vaf = as.numeric(s$vaf),
        source_id = paste0(s$chrom, ":", s$pos, ":", s$ref, ">", s$alt), genes = toupper(s$gene),
        dp = NA_real_, callers = "", locus = paste0(s$chrom, ":", kf_int(s$pos)),
        loh = FALSE, in_panel = TRUE)
    }
  }

  cn = hmf$cn
  if (!is.null(cn) && nrow(cn) > 0) {
    c2 = copy(cn)
    if (!"sv_size" %in% names(c2)) c2[, sv_size := NA_real_]
    sv_src = grepl(" SV$", c2$source)
    c2[, tier := fifelse(!sv_src, fifelse(call == "PARTIAL_AMP", 2L, 1L),
                         fifelse(!is.na(sv_size) & sv_size <= KF_FOCAL_SV_MAX, 2L, 3L))]
    c2[, dir := fifelse(call == "DEL", "DEL", "AMP")]
    g = c2[, {
      best = which.min(tier)
      cnv = c(min_cn, max_cn); cnv = cnv[!is.na(cnv)]
      cn_txt = if (length(cnv) == 0) "" else
               if (min(cnv) == max(cnv)) trimws(formatC(min(cnv), format = "fg", digits = 3))
               else paste0(trimws(formatC(min(cnv), format = "fg", digits = 3)), "–",
                           trimws(formatC(max(cnv), format = "fg", digits = 3)))
      svn = note[grepl(" SV$", source) & nzchar(note)]
      .(call = call[best], tier = min(tier), chrom = chrom[1], start = start[1], end = end[1],
        role = role[1], sources = paste(unique(source), collapse = ", "), cn_txt = cn_txt,
        sv_note = if (length(svn)) svn[1] else "")
    }, by = .(gene, dir)]
    rows$cn = data.table(
      gene = toupper(g$gene), type = "CNA",
      alteration = paste0(g$call, ifelse(nzchar(g$cn_txt), paste0(" (CN ", g$cn_txt, ")"),
                                         " (SV deletion)")),
      evidence = kf_join(paste0("HMF ", g$call, " · ", g$sources), g$sv_note,
                         ifelse(g$tier == 3L, "large SV deletion, possibly heterozygous", "")),
      tier = as.integer(g$tier), role = na0(g$role), panels = tag, vaf = NA_real_,
      source_id = paste0(g$chrom, ":", g$start, "-", g$end), genes = toupper(g$gene),
      dp = NA_real_, callers = g$sources,
      locus = paste0(g$chrom, ":", kf_int(g$start), "–", kf_int(g$end)),
      loh = FALSE, in_panel = TRUE)
  }

  fu = hmf$fusions
  if (!is.null(fu) && nrow(fu) > 0) {
    f = fu[reported %in% TRUE]
    if (nrow(f) > 0) {
      f[, caller_lab := fifelse(caller %in% names(KF_CALLER_LABEL), KF_CALLER_LABEL[caller], na0(caller))]
      g = f[, .(five_exon = five_exon[1], three_exon = three_exon[1], frame = frame[1],
                likelihood = min(likelihood), high_impact = any(high_impact %in% TRUE),
                callers = paste(unique(caller_lab[nzchar(caller_lab)]), collapse = ", "),
                note = { n = note[nzchar(na0(note))]; if (length(n)) n[1] else "" },
                ids = paste(unique(sv_id), collapse = ","), locus = locus[1]),
            by = .(type, five_gene, three_gene)]
      members = lapply(seq_len(nrow(g)), function(i)
        unique(toupper(na.omit(c(g$five_gene[i], g$three_gene[i])))))
      exons = kf_join(ifelse(is.na(g$five_exon), "", paste0("5′ ", g$five_exon)),
                      ifelse(is.na(g$three_exon), "", paste0("3′ ", g$three_exon)), sep = ", ")
      rows$fusions = data.table(
        gene = vapply(members, paste, character(1), collapse = "–"), type = "SV",
        alteration = paste0(gsub("_", " ", tolower(g$type)), " fusion",
                            ifelse(nzchar(exons), paste0(" (", exons, ")"), "")),
        evidence = kf_join(paste0("HMF ", g$type), ifelse(g$frame == "n/a", "", g$frame),
                           ifelse(g$high_impact, "high-impact partner", ""), g$note,
                           ifelse(nzchar(g$callers), paste0("called by ", g$callers), "")),
        tier = ifelse(g$likelihood == "HIGH", 1L, 2L), role = "", panels = tag, vaf = NA_real_,
        source_id = g$ids, genes = vapply(members, paste, character(1), collapse = ","),
        dp = NA_real_, callers = g$callers, locus = g$locus, loh = FALSE, in_panel = TRUE)
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

# ---- findings ------------------------------------------------------------------------------

build_findings = function(variant_table, section_data, all_panels, reference, hmf = NULL) {
  panels = if (is.null(all_panels)) list() else all_panels[vapply(all_panels, is.list, logical(1))]
  # HMF-built lists take their findings from build_hmf() (kf_hmf_rows()); the generic rules run
  # over the other lists. Without an HMF result (it failed to load) the list falls back to them.
  hmf_lists = names(panels)[vapply(panels, kf_is_hmf_panel, logical(1))]
  if (is.null(hmf) || !is.null(hmf$error)) hmf_lists = character(0)
  all_list_names = names(panels)
  panels = panels[setdiff(names(panels), hmf_lists)]
  roles = kf_gene_roles(panels)
  snv_sets = lapply(panels, kf_panel_snv)
  cna = section_data$cna
  loh_genes = if (!is.null(cna) && !is.null(cna$genes) && nrow(cna$genes) > 0 &&
                  "loh" %in% names(cna$genes))
                unique(toupper(cna$genes$gene[cna$genes$loh %in% TRUE])) else character(0)

  sv = kf_small_variant_rows(variant_table, snv_sets, roles, loh_genes)
  n_nonpanel = attr(sv, "n_nonpanel")
  sv_tab = if (!is.null(section_data$sv)) section_data$sv$table else NULL
  cna_rows = kf_cna_rows(cna, roles, reference)
  if (length(hmf_lists) && nrow(cna_rows)) {
    # The copy-number section's gene calls carry every list's name, the HMF one included
    cna_rows[, panels := vapply(strsplit(panels, ",", fixed = TRUE),
                                function(p) paste(setdiff(p, hmf_lists), collapse = ","), character(1))]
    cna_rows = cna_rows[nzchar(panels)]
  }
  rows = rbindlist(list(empty_findings(), sv,
                        kf_sv_rows(sv_tab, panels, roles),
                        cna_rows,
                        kf_hmf_rows(hmf, hmf_lists)), use.names = TRUE, fill = TRUE)
  if (nrow(rows)) {
    rows[is.na(rescued), rescued := FALSE]
    rows[is.na(ev_tier), ev_tier := tier]
  }
  # LOH annotates every row of its gene
  if (nrow(rows)) rows[gene %in% loh_genes, loh := TRUE]
  # One group per gene (or per joined rearrangement): its rows are contiguous, groups ordered
  # by their best tier, rows within a group by tier and type
  rows[, type_rank := match(type, KF_TYPES)]
  rows[, group_tier := tier]
  if (nrow(rows)) rows[, group_tier := min(tier), by = gene]
  setorder(rows, group_tier, gene, tier, type_rank, -vaf, na.last = TRUE)
  rows[, c("type_rank", "group_tier") := NULL]
  per_panel_n = vapply(all_list_names, function(nm) {
    m = rows$tier <= 2L & vapply(strsplit(rows$panels, ",", fixed = TRUE),
                                 function(p) nm %in% p, logical(1))
    length(unique(rows$gene[m]))
  }, integer(1))
  list(findings = rows[], per_panel_n = per_panel_n, n_nonpanel = n_nonpanel)
}

# ---- circos driver labels ------------------------------------------------------------------

# Tiers the genome circos labels (assets/js/circos_drivers.js shows those of the ticked lists)
CIRCOS_LABEL_TIERS = 1:2

# Tier 1-2 finding lines as label rows: one row per (line, member gene), so the joined
# "EWSR1–ATF1" row labels EWSR1 and ATF1 at their own loci. `group` is the line's `gene`
# (its key-findings group), `genes` the group's members, which is what a Custom… list matches.
# A whitelist-rescued line whose evidence is tier 1-2 is a row too, at that tier with
# `rescued` 1: the browser labels it only while rescued SVs are shown.
circos_driver_rows = function(findings) {
  empty = data.table(label = character(), group = character(), genes = character(),
                     panels = character(), tier = integer(), type = character(),
                     alteration = character(), rescued = integer())
  if (is.null(findings) || nrow(findings) == 0) return(empty)
  resc = if ("rescued" %in% names(findings)) findings$rescued %in% TRUE else rep(FALSE, nrow(findings))
  evt = if ("ev_tier" %in% names(findings)) findings$ev_tier else findings$tier
  # Non-panel rows (panels "") stay: a typed Custom… list labels them, as the overview shows them
  keep = (!resc & findings$tier %in% CIRCOS_LABEL_TIERS) | (resc & evt %in% CIRCOS_LABEL_TIERS)
  f = findings[keep]
  if (nrow(f) == 0) return(empty)
  f_tier = ifelse(resc[keep], evt[keep], f$tier)
  members = strsplit(ifelse(nzchar(f$genes), f$genes, f$gene), ",", fixed = TRUE)
  idx = rep(seq_len(nrow(f)), lengths(members))
  out = data.table(label = trimws(unlist(members)), group = f$gene[idx], genes = f$genes[idx],
                   panels = f$panels[idx], tier = as.integer(f_tier[idx]), type = f$type[idx],
                   alteration = f$alteration[idx], rescued = as.integer(resc[keep][idx]))
  unique(out[nzchar(label)])
}

# One locus per label gene, data.table(gene, chrom, pos): the midpoint of its first interval in
# any embedded list (every row with coordinates, not only SV-scoped ones), else its gene-level
# CN row, else its first small variant in the table; a gene with none of these has no label.
circos_label_loci = function(genes, all_panels, section_data = list(), variant_table = NULL) {
  empty = data.table(gene = character(), chrom = character(), pos = numeric())
  genes = unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) == 0) return(empty)
  key = toupper(genes)
  panels = if (is.null(all_panels)) list() else all_panels[vapply(all_panels, is.list, logical(1))]
  src = list(rbindlist(lapply(panels, panel_gene_intervals, scope = "all")))
  cg = section_data$cna$genes
  if (!is.null(cg) && nrow(cg) > 0 && all(c("gene", "chrom", "start", "end") %in% names(cg)))
    src[[2]] = cg[, .(chrom, start, end, gene)]
  if (!is.null(variant_table) && nrow(variant_table) > 0 && "symbol" %in% names(variant_table))
    src[[3]] = variant_table[!is.na(symbol), .(chrom, start = pos, end = pos, gene = symbol)]
  src = rbindlist(src, use.names = TRUE, fill = TRUE)
  if (nrow(src) == 0) return(empty)
  src = src[!is.na(chrom) & !is.na(start) & !is.na(end)]
  src[, gkey := toupper(gene)]
  src = src[gkey %in% key & !duplicated(gkey)]
  if (nrow(src) == 0) return(empty)
  data.table(gene = genes[match(src$gkey, key)], chrom = ensure_chr_prefix(src$chrom),
             pos = (as.numeric(src$start) + as.numeric(src$end)) / 2)
}

# ---- biomarkers ----------------------------------------------------------------------------

kf_bm = function(id, label, value, detail = NULL, status = "ok") {
  list(id = id, label = label, value = value,
       detail = if (length(detail) && any(nzchar(detail))) paste(detail[nzchar(detail)], collapse = " · ") else "",
       status = status)
}

kf_num = function(x, digits = 2) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  formatC(as.numeric(x), format = "f", digits = digits)
}

kf_pct = function(x, digits = 1) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  paste0(formatC(100 * as.numeric(x), format = "f", digits = digits), "%")
}

# Purity disagreement between the callers' top solutions (README "Copy-number calls")
PURITY_SPREAD_WARN      = 0.2   # --purity-spread-warn default: max - min purity
PURITY_WGD_PLOIDY_RATIO = 1.6   # max / min ploidy at which the callers disagree about WGD
PURITY_SOURCE_LABELS    = c(wakhan = "Wakhan", ascat = "ASCAT", savana = "SAVANA")

# `pbs` is cna$purity_by_source, already one top solution per caller (Wakhan rank 1, ASCAT's
# fit, SAVANA rank 1). A caller without a purity drops out; NULL with fewer than two left.
purity_disagreement = function(pbs, threshold = PURITY_SPREAD_WARN, used = NA_character_) {
  if (is.null(pbs) || !is.data.frame(pbs) || nrow(pbs) == 0 ||
      !all(c("source", "purity") %in% names(pbs))) return(NULL)
  if (is.null(threshold) || length(threshold) == 0 || is.na(threshold)) threshold = PURITY_SPREAD_WARN
  num = function(x) suppressWarnings(as.numeric(x))
  d = data.table(source = as.character(pbs$source), purity = num(pbs$purity),
                 ploidy = if ("ploidy" %in% names(pbs)) num(pbs$ploidy) else NA_real_)
  d = d[!is.na(purity) & !duplicated(source)]
  if (nrow(d) < 2) return(NULL)
  used = if (length(used) == 1 && !is.na(used) && used %in% d$source) used else NA_character_
  if (!is.na(used)) d = d[order(source != used)]
  spread = round(max(d$purity) - min(d$purity), 2)
  pl = d$ploidy[!is.na(d$ploidy) & d$ploidy > 0]
  ratio = if (length(pl) > 1) max(pl) / min(pl) else NA_real_
  list(sources = d, spread = spread, threshold = as.numeric(threshold),
       flag = spread >= threshold - 1e-9,
       ploidy_range = if (length(pl) > 1) range(pl) else c(NA_real_, NA_real_),
       ploidy_ratio = ratio,
       wgd_mismatch = !is.na(ratio) && ratio >= PURITY_WGD_PLOIDY_RATIO,
       used = used)
}

# The source whose purity/ploidy the CN calls use: the segment source, else the `used` row
purity_used_source = function(cna) {
  if (!is.null(cna$source) && length(cna$source) == 1 && !is.na(cna$source)) return(cna$source)
  pbs = cna$purity_by_source
  if (is.null(pbs) || nrow(pbs) == 0 || !"used" %in% names(pbs)) return(NA_character_)
  u = pbs$source[pbs$used %in% TRUE]
  if (length(u)) u[1] else NA_character_
}

# The callout naming every caller's purity and ploidy; NULL unless `pc$flag`
purity_disagreement_notice = function(pc, verdict = NULL, tumour_only = FALSE) {
  if (is.null(pc) || !isTRUE(pc$flag)) return(NULL)
  lab = function(x) ifelse(x %in% names(PURITY_SOURCE_LABELS), PURITY_SOURCE_LABELS[x], x)
  s = pc$sources
  each = paste0(lab(s$source), " ", vapply(s$purity, kf_num, ""), " (ploidy ",
                vapply(s$ploidy, kf_num, ""), ")")
  used = if (!is.na(pc$used)) s[source == pc$used] else NULL
  v_src = if (!is.null(verdict$purity_source) && !is.na(verdict$purity_source)) verdict$purity_source
          else NA_character_
  v_pur = if (!is.null(verdict$purity)) verdict$purity else NA_real_
  v_known = !is.na(v_src) && v_src %in% names(PURITY_SOURCE_LABELS) && !is.na(v_pur)
  verdict_txt = if (!isTRUE(tumour_only) || identical(verdict$status, "not_run")) NULL
    else if (identical(verdict$status, "skipped"))
      paste0(" ClairS-TO germline tagging (Verdict) was skipped for this sample",
             if (v_known) paste0(" (", lab(v_src), "'s purity, ", kf_num(v_pur), ")") else "", ".")
    else if (identical(verdict$status, "applied") && v_known)
      paste0(" ClairS-TO germline tagging (Verdict) used ", lab(v_src), "'s purity, ", kf_num(v_pur), ".")
    else " In a tumour-only run, ClairS-TO germline tagging (Verdict) depends on purity too."
  tags$div(class = "purity-alert", role = "note",
    tags$p(class = "purity-alert__title",
           tags$span(class = "purity-alert__icon", `aria-hidden` = "true", "⚠"),
           paste0("Purity estimates disagree by ", kf_num(pc$spread))),
    tags$p(class = "purity-alert__values", paste(each, collapse = " · "),
           paste0(" — top solution of each; this warning shows at a difference of ",
                  kf_num(pc$threshold), " or more.")),
    if (isTRUE(pc$wgd_mismatch))
      tags$p(paste0("Their ploidies differ ", formatC(pc$ploidy_ratio, format = "f", digits = 1),
                    "-fold (", kf_num(pc$ploidy_range[1]), " to ", kf_num(pc$ploidy_range[2]),
                    "), so the callers also disagree about whole-genome doubling.")),
    tags$p(paste0("Copy-number calls and the copy-number key findings depend on purity and ploidy",
                  if (!is.null(used) && nrow(used))
                    paste0("; this report uses ", lab(used$source), "'s (", kf_num(used$purity),
                           " / ", kf_num(used$ploidy), ").") else ".",
                  " Check the copy-number profile before relying on them.",
                  if (!is.null(verdict_txt)) verdict_txt else "")))
}

build_biomarkers = function(section_data, ascat_pp, tmb_info, counts, purity_check = NULL) {
  src_lab = c(wakhan = "Wakhan", ascat = "ASCAT", savana = "SAVANA")
  lab = function(x) { l = src_lab[x]; ifelse(is.na(l), x, l) }
  cna = section_data$cna
  out = list()

  # Purity & ploidy: the source the CN calls use, named; the others small
  pur = NA_real_; plo = NA_real_; src = NA_character_; others = character(0)
  if (!is.null(cna) && (!is.null(cna$purity) || !is.null(cna$ploidy))) {
    pur = as.numeric(cna$purity %||kf% NA); plo = as.numeric(cna$ploidy %||kf% NA)
    pbs = cna$purity_by_source
    src = purity_used_source(cna)
    if (!is.null(pbs) && nrow(pbs) > 0) {
      o = pbs[!(source %in% src)]
      others = sprintf("%s %s / %s", lab(o$source), vapply(o$purity, kf_num, ""),
                       vapply(o$ploidy, kf_num, ""))
      if (is.null(purity_check)) purity_check = purity_disagreement(pbs, used = src)
    }
  }
  if (is.na(pur) && is.na(plo) && !is.null(ascat_pp)) {
    pur = as.numeric(ascat_pp$purity %||kf% NA); plo = as.numeric(ascat_pp$ploidy %||kf% NA)
    if (!is.na(pur) || !is.na(plo)) src = "ascat"
  }
  if (is.na(pur) && is.na(plo)) {
    out$purity = kf_bm("purity", "Purity / ploidy", "not run",
                       "no ASCAT, Wakhan or SAVANA estimate found", "not_run")
  } else {
    disagree = isTRUE(purity_check$flag)
    out$purity = kf_bm("purity", "Purity / ploidy",
                       paste0(kf_num(pur), " / ", kf_num(plo)),
                       c(if (!is.na(src)) lab(src) else "",
                         if (length(others)) paste0("also ", paste(others, collapse = "; ")) else "",
                         if (disagree) paste0("purity estimates differ by ", kf_num(purity_check$spread))
                         else ""),
                       if (disagree) "warn" else "ok")
  }

  # Whole-genome duplication
  gen = if (!is.null(cna)) cna$genome else NULL
  if (is.null(gen) || is.null(gen$wgd) || is.na(gen$wgd)) {
    aw = if (!is.null(gen)) gen$ascat_wgd else NULL
    if (!is.null(aw) && length(aw) && !is.na(aw)) {
      out$wgd = kf_bm("wgd", "Whole-genome duplication", if (aw > 0) "yes" else "no", "ASCAT")
    } else out$wgd = kf_bm("wgd", "Whole-genome duplication", "not run",
                           "no copy-number segments", "not_run")
  } else {
    basis = c(
      if (!is.null(gen$major2_fraction) && !is.na(gen$major2_fraction))
        sprintf("%s: major allele ≥ 2 over %s of autosomes", lab(gen$wgd_source %||kf% ""),
                kf_pct(gen$major2_fraction, 0)) else lab(gen$wgd_source %||kf% ""),
      if (!is.null(gen$ascat_wgd) && !is.na(gen$ascat_wgd) &&
          (gen$ascat_wgd > 0) != isTRUE(gen$wgd))
        paste0("ASCAT says ", if (gen$ascat_wgd > 0) "yes" else "no") else "")
    out$wgd = kf_bm("wgd", "Whole-genome duplication", if (isTRUE(gen$wgd)) "yes" else "no",
                    basis, if (any(grepl("ASCAT says", basis))) "warn" else "ok")
  }

  # LOH fraction
  if (!is.null(gen) && !is.null(gen$loh_fraction) && !is.na(gen$loh_fraction)) {
    out$loh = kf_bm("loh", "LOH fraction", kf_pct(gen$loh_fraction),
                    c(paste0("of autosomes, ", lab(gen$wgd_source %||kf% "")),
                      if (!is.null(gen$ascat_loh) && !is.na(gen$ascat_loh) && !identical(gen$wgd_source, "ascat"))
                        paste0("ASCAT ", kf_pct(gen$ascat_loh)) else ""))
  } else out$loh = kf_bm("loh", "LOH fraction", "not run", "no copy-number segments", "not_run")

  # Coding TMB (compute_tmb()): SNVs lead, since long-read homopolymer indel errors can make up
  # most of the total; the total and the split follow, all over unique variants
  if (!is.null(tmb_info) && !is.null(tmb_info$tmb) && !is.na(tmb_info$tmb)) {
    has_split = !is.null(tmb_info$tmb_snv) && !is.na(tmb_info$tmb_snv)
    out$tmb = kf_bm("tmb", "Coding TMB",
                    if (has_split) paste0(kf_num(tmb_info$tmb_snv), " mut/Mb (SNVs)")
                    else paste0(kf_num(tmb_info$tmb), " mut/Mb"),
                    if (!is.null(tmb_info$n_nonsyn) && !is.na(tmb_info$n_nonsyn))
                      if (has_split)
                        paste0(kf_int(tmb_info$n_snv), " SNV", if (tmb_info$n_snv == 1) "" else "s",
                               " + ", kf_int(tmb_info$n_indel), " indel",
                               if (tmb_info$n_indel == 1) "" else "s", " non-syn / ",
                               tmb_info$denominator_mb, " Mb · all ", kf_num(tmb_info$tmb), " mut/Mb")
                      else paste0(kf_int(tmb_info$n_nonsyn), " non-syn / ", tmb_info$denominator_mb, " Mb"))
  } else out$tmb = kf_bm("tmb", "Coding TMB", "not run", "no small-variant calls", "not_run")

  # Small variants
  smv = counts$small_variants
  if (!is.null(smv) && !is.null(smv$total) && !is.na(smv$total)) {
    out$small = kf_bm("small", "Small variants", kf_int(smv$total),
                      paste0(kf_int(smv$snv), " SNVs · ", kf_int(smv$indel), " indels"))
  } else out$small = kf_bm("small", "Small variants", "not run", "no VEP output", "not_run")

  # SVs per caller: never summed, two callers report the same rearrangements
  svd = section_data$sv
  nbc = if (!is.null(svd) && !is.null(svd$n_by_caller)) svd$n_by_caller
        else if (!is.null(svd) && !is.null(svd$table) && nrow(svd$table) > 0 && "caller" %in% names(svd$table))
          { tb = table(svd$table$caller); setNames(as.integer(tb), names(tb)) }
        else NULL
  if (length(nbc) > 0) {
    cl = KF_CALLER_LABEL[names(nbc)]; cl[is.na(cl)] = names(nbc)[is.na(cl)]
    # Whitelist-rescued rows are counted apart: "Severus 323 (+9,388 whitelist-rescued)"
    nres = svd$n_rescued_by_caller
    nres = if (is.null(nres)) setNames(rep(0L, length(nbc)), names(nbc))
           else setNames(as.integer(nres[names(nbc)]), names(nbc))
    nres[is.na(nres)] = 0L
    tbl = svd$table
    pass_rows = if (!is.null(tbl) && nrow(tbl) > 0) !sv_is_rescued(tbl) else logical(0)
    tr = if (!is.null(tbl) && nrow(tbl) > 0 && all(c("caller", "svclass") %in% names(tbl)))
           vapply(names(nbc), function(c) sum(pass_rows & tbl$caller == c & tbl$svclass == "translocation"), integer(1))
         else NULL
    out$sv = kf_bm("sv", "Structural variants",
                   paste(sprintf("%s %s%s", cl, kf_int(nbc - nres),
                                 ifelse(nres > 0, paste0(" (+", kf_int(nres), " whitelist-rescued)"), "")),
                         collapse = " · "),
                   c(if (!is.null(tr)) paste0("translocations: ", paste(sprintf("%s %s", cl, kf_int(tr)), collapse = ", ")) else "",
                     if (length(nbc) > 1) "per caller, not summed" else ""))
  } else out$sv = kf_bm("sv", "Structural variants", "not run", "no SV caller output", "not_run")

  # Chromothripsis: one card per ShatterSeek minimum-tier setting (the report's selector)
  ss = section_data$shatterseek
  if (is.null(ss) || !identical(ss$status, "ok")) {
    out$chromothripsis = kf_bm("chromothripsis", "Chromothripsis", "not run", "ShatterSeek", "not_run")
  } else {
    tiers = as.integer(names(SHATTERSEEK_MIN_TIERS))
    views = lapply(tiers, function(t) kf_chromothripsis_bm(shatterseek_at_tier(ss, t)))
    names(views) = tiers
    out$chromothripsis = views[[as.character(SHATTERSEEK_MIN_TIER_DEFAULT)]]
    if (!is.null(ss$high) && nrow(ss$high) > 0) out$chromothripsis$views = views
  }

  # ecDNA / focal amplicons (detected != ecDNA: read n_ecdna/n_bfb)
  ec = section_data$ecdna
  st = if (is.null(ec) || is.null(ec$status)) "not_run" else ec$status
  n_amp = if (!is.null(ec$amplicons)) nrow(ec$amplicons) else 0L
  out$ecdna = switch(st,
    not_run = kf_bm("ecdna", "ecDNA / amplicons", "not run", "CoRAL / AmpliconClassifier", "not_run"),
    none_detected = kf_bm("ecdna", "ecDNA / amplicons", "none detected", ecdna_verdict(ec)),
    unclassified = kf_bm("ecdna", "ecDNA / amplicons", "focal amplicons found, not classified",
                         ecdna_verdict(ec), "warn"),
    detected = kf_bm("ecdna", "ecDNA / amplicons",
                     { bits = c(if (isTRUE(ec$n_ecdna > 0)) sprintf("%d ecDNA+", ec$n_ecdna),
                                if (isTRUE(ec$n_bfb > 0)) sprintf("%d BFB+", ec$n_bfb))
                       if (length(bits)) paste(bits, collapse = ", ")
                       else "no ecDNA / BFB" },
                     ecdna_verdict(ec)),
    kf_bm("ecdna", "ecDNA / amplicons", "unknown", ecdna_verdict(ec), "warn"))

  # Dominant signatures: only `dominant` (well-fitted contexts, >= SIG_DOMINANT_MIN_PCT)
  sg = section_data$signatures
  if (is.null(sg)) {
    out$signatures = kf_bm("signatures", "Dominant signatures", "not run",
                           "SigProfilerAssignment", "not_run")
  } else if (is.null(sg$dominant) || nrow(sg$dominant) == 0 ||
             !any(sg$dominant$context == KF_SIGNATURE_CONTEXT)) {
    out$signatures = kf_bm("signatures", "Dominant signatures", "no well-fitted signature",
                           "no SBS fit passed the fit / mutation-count thresholds")
  } else {
    # SBS only: indel (ID) and doublet (DBS) fits stay in the signatures section. On long
    # reads an ID fit can pass the cosine threshold while resting on homopolymer indel
    # errors, and "Defective MMR" on the first screen from that is a claim the data cannot carry.
    d = sg$dominant[sg$dominant$context == KF_SIGNATURE_CONTEXT]
    art = if ("artefact" %in% names(d)) d$artefact %in% TRUE else rep(FALSE, nrow(d))
    aet = if ("aetiology" %in% names(d)) as.character(d$aetiology) else rep(NA_character_, nrow(d))
    items = paste0(d$signature, " ", formatC(as.numeric(d$pct), format = "f", digits = 0), "%",
                   ifelse(art, " [artefact]", ""),
                   ifelse(is.na(aet) | !nzchar(aet), "", paste0(" (", aet, ")")))
    out$signatures = kf_bm("signatures", "Dominant signatures",
                           paste(d$signature[!art][seq_len(min(3, sum(!art)))], collapse = ", "),
                           items)
    if (!nzchar(out$signatures$value)) out$signatures$value = "artefact only"
    out$signatures$items = items
  }
  out
}

# The chromothripsis card for one ShatterSeek result: the called chromosomes, each with the PCAWG
# tier(s) it meets and the caller combinations that called it
kf_chromothripsis_bm = function(ss) {
  if (is.null(ss$high) || nrow(ss$high) == 0) {
    below = if (isTRUE(ss$n_below_tier > 0))
              sprintf("%d call%s below this tier", ss$n_below_tier, if (ss$n_below_tier == 1) "" else "s") else ""
    return(kf_bm("chromothripsis", "Chromothripsis", if (nzchar(below)) "none at this tier" else "none called",
                 c(paste0("ShatterSeek: ", paste(unname(ss$labels), collapse = ", ")), below)))
  }
  h = ss$high
  by = h[, .(combos = paste(unique(combo_label), collapse = ", "),
             tiers = shatterseek_tiers(criteria),
             conf = paste(unique(confidence), collapse = "/")), by = chrom]
  dis = isTRUE(shatterseek_disagree(ss))
  det = paste0(by$chrom, " ", ifelse(nzchar(by$tiers), by$tiers, by$conf), " (", by$combos, ")")
  kf_bm("chromothripsis", "Chromothripsis", paste(by$chrom, collapse = ", "),
        c(det, if (dis) paste0("caller combinations disagree; ran: ",
                                paste(unname(ss$labels), collapse = ", ")) else ""),
        if (dis) "warn" else "ok")
}

# Signature context the overview reports; the others are shown in their own section only
KF_SIGNATURE_CONTEXT = "SBS"

`%||kf%` = function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- entry point ---------------------------------------------------------------------------

build_overview = function(variant_table, section_data, all_panels, ascat_pp, tmb_info,
                          counts, ctx, hmf = NULL) {
  if (is.null(section_data)) section_data = list()
  if (is.null(counts)) counts = list()
  if (is.null(ctx)) ctx = list()
  ref = if (!is.null(ctx$reference)) ctx$reference else NA_character_
  f = build_findings(variant_table, section_data, all_panels, ref, hmf)
  mode = if (!is.null(ctx$outputs) && !is.null(ctx$outputs$mode)) ctx$outputs$mode else NA_character_
  cna = section_data$cna
  pc = if (!is.null(cna)) purity_disagreement(cna$purity_by_source, ctx$purity_spread_warn,
                                              purity_used_source(cna))
  list(findings    = f$findings,
       biomarkers  = build_biomarkers(section_data, ascat_pp, tmb_info, counts, pc),
       purity_check = pc,
       per_panel_n = f$per_panel_n,
       tumour_only = identical(as.character(mode), "tumour-only"),
       verdict     = if (!is.null(ctx$outputs)) parse_verdict_status(ctx$outputs$verdict_status) else NULL,
       n_nonpanel  = f$n_nonpanel)
}
