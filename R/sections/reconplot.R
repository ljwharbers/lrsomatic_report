# ReConPlot figures section module.
#
# Inputs: reconplot/<combo>/genome_wide/<id>_genome_wide.png,
#         reconplot/<combo>/per_chromosome/<id>_chr*.png and reconplot/<combo>/focus/<id>_*.png
#         (the user's --reconplot_regions / --reconplot_genes), combo in severus_ascat,
#         severus_wakhan, savana (any other combo directory is kept, after those). PDFs and the
#         harmonised TSVs are ignored.
#
# Contract
#   id "reconplot"; parse(inputs, section_data, ctx) returns NULL (no figure at all) or
#   list(combos = list(<combo> = list(genome = <path or NULL>,
#                                     chroms = c(chr1 = <path>, ...),     # named chr vector
#                                     focus  = c(<region label> = <path>, ...))))
#   chroms is in display order (ctx$chromosomes, unknown contigs last); focus in file-name
#   order, labelled by the file stem without the sample prefix. Rendered by
#   templates/sections/_reconplot.qmd (included at the end of _sv.qmd).

# Which per-chromosome figures slim mode embeds, and why. One row per chromosome:
#   data.table(chrom, always, panels, genes, reasons, kinds, ss_tier, priority)
#   always   TRUE for a list-independent reason: a ShatterSeek call (any combination, any
#            confidence but "none") or an SV with a Padfoot cancer-gene rating: any HIGH_ class
#            (HIGH_ONCO, HIGH_TSG; the set is open-ended, as _sv.qmd also treats it), on either
#            breakend's contig, unless the SV was whitelist-rescued (sv_is_rescued()).
#            Plain HIGH (score 3) does not count: on CCS15-PB it
#            selected 12 chromosomes for USH2A, PRSS2, MACROD2 and other fragile-site genes;
#            HIGH_ONCO picked out exactly ALK, CDKN2A, EWSR1–ATF1, NSD3, PRKN and CNTNAP2
#   panels / genes  comma-joined lists and member genes of the tier 1-2 SV / CNA key findings
#            on the chromosome (OVERVIEW$findings); the browser shows such a figure only while
#            one of those lists is ticked, or a Custom… list names one of those genes
#   reasons  short caption text, e.g. "EWSR1–ATF1 translocation (tier 1)"
#   kinds    comma-joined reason kinds: chromothripsis, padfoot, finding (the report's selector)
#   ss_tier  best PCAWG tier of the chromosome's chromothripsis calls (NA: none, or none named)
#   priority 1 chromothripsis, 2 Padfoot, 3 tier-1 finding, 4 tier-2 finding: the best reason,
#            the order in which figures are kept under the slim-mode cap
# The genome-wide figure is not selectable: at 10200x1200 px it is unreadable at page width.
RECONPLOT_FINDING_TYPES = c("SV", "CNA")
RECONPLOT_PADFOOT_IMPACT = "^HIGH_"
RECONPLOT_MAX_REASONS   = 3L

reconplot_chromosomes = function(findings = NULL, sv_table = NULL, shatterseek_high = NULL) {
  out = list()
  if (!is.null(shatterseek_high) && nrow(shatterseek_high) > 0 && "chrom" %in% names(shatterseek_high)) {
    keep = if ("confidence" %in% names(shatterseek_high))
             !is.na(shatterseek_high$confidence) & tolower(shatterseek_high$confidence) != "none"
           else rep(TRUE, nrow(shatterseek_high))
    ssd = data.table(chrom = ensure_chr_prefix(as.character(shatterseek_high$chrom[keep])),
                     criteria = if ("criteria" %in% names(shatterseek_high))
                                  as.character(shatterseek_high$criteria[keep]) else NA_character_)
    if (nrow(ssd)) {
      ssd = ssd[, .(tiers = shatterseek_tiers(criteria),
                    ss_tier = suppressWarnings(min(shatterseek_best_tier(criteria), na.rm = TRUE))), by = chrom]
      ssd[!is.finite(ss_tier), ss_tier := NA_real_]
      out[[length(out) + 1]] = data.table(
        chrom = ssd$chrom, always = TRUE, panels = "", genes = "", kind = "chromothripsis",
        ss_tier = as.integer(ssd$ss_tier), prio = 1L,
        reason = paste0("chromothripsis call", ifelse(nzchar(ssd$tiers), paste0(" (", ssd$tiers, ")"), "")))
    }
  }
  if (!is.null(sv_table) && nrow(sv_table) > 0 && "padfoot_impact" %in% names(sv_table)) {
    # A whitelist-rescued call is hidden by default everywhere else, so it selects nothing here
    hi = sv_table[!sv_is_rescued(sv_table) & !is.na(padfoot_impact) &
                  grepl(RECONPLOT_PADFOOT_IMPACT, padfoot_impact)]
    if (nrow(hi) > 0) {
      # Padfoot genes read "ATF1 (intron3), EWSR1 (intron7)": symbols only, one reason per
      # rating, so mates and repeated calls in one gene collapse
      g = if ("padfoot_genes" %in% names(hi)) as.character(hi$padfoot_genes) else rep("", nrow(hi))
      g[is.na(g)] = ""
      g = gsub("\\s*\\([^)]*\\)", "", g)
      pf = rbind(data.table(chrom = as.character(hi$chrom_a), impact = hi$padfoot_impact, genes = g),
                 data.table(chrom = as.character(hi$chrom_b), impact = hi$padfoot_impact, genes = g))
      pf = pf[!is.na(chrom) & nzchar(chrom)]
      if (nrow(pf)) {
        pf = pf[, .(genes = {
                      s = unique(trimws(unlist(strsplit(genes, ",", fixed = TRUE))))
                      paste(s[nzchar(s)], collapse = ", ")
                    }), by = .(chrom = ensure_chr_prefix(chrom), impact)]
        out[[length(out) + 1]] = data.table(
          chrom = pf$chrom, always = TRUE, panels = "", genes = "", kind = "padfoot",
          ss_tier = NA_integer_, prio = 2L,
          reason = paste0("Padfoot ", pf$impact, ifelse(nzchar(pf$genes), paste0(": ", pf$genes), "")))
      }
    }
  }
  if (!is.null(findings) && nrow(findings) > 0) {
    f = findings[tier <= 2L & type %in% RECONPLOT_FINDING_TYPES]
    if (nrow(f) > 0) {
      # Loci are "chr12–chr22", "chr9:21.90–22.01 Mb" or "chr9:21,967,752–21,995,301"
      chs = regmatches(f$locus, gregexpr("chr[A-Za-z0-9_.]+", f$locus))
      idx = rep(seq_len(nrow(f)), lengths(chs))
      if (length(idx)) {
        # The alteration up to its locus or caller detail: "EWSR1–ATF1 translocation",
        # "CDKN2A deletion", "homozygous deletion" (CN lines do not name their gene)
        alt = trimws(sub("\\s+(chr[A-Za-z0-9_.]+[:–]|\\().*$", "", as.character(f$alteration[idx])))
        gene = f$gene[idx]
        what = ifelse(startsWith(alt, gene), alt, paste(gene, alt))
        out[[length(out) + 1]] = data.table(
          chrom = unlist(chs), always = FALSE, panels = f$panels[idx],
          genes = ifelse(nzchar(f$genes[idx]), f$genes[idx], f$gene[idx]), kind = "finding",
          ss_tier = NA_integer_, prio = 2L + as.integer(f$tier[idx]),
          reason = paste0(what, " (tier ", f$tier[idx], ")"))
      }
    }
  }
  empty = data.table(chrom = character(), always = logical(), panels = character(),
                     genes = character(), reasons = character(), kinds = character(),
                     ss_tier = integer(), priority = integer())
  if (length(out) == 0) return(empty)
  d = rbindlist(out)
  join = function(x) { x = unique(unlist(strsplit(x, ",", fixed = TRUE))); paste(x[nzchar(x)], collapse = ",") }
  # The lists of the finding reasons are kept on an always-shown chromosome too: the reason
  # selector can switch its list-independent reasons off
  d[, .(always  = any(always),
        panels  = join(panels),
        genes   = join(genes),
        kinds   = paste(unique(kind), collapse = ","),
        ss_tier = if (any(!is.na(ss_tier))) min(ss_tier, na.rm = TRUE) else NA_integer_,
        priority = min(prio),
        reasons = {
          r = unique(reason)
          if (length(r) > RECONPLOT_MAX_REASONS)
            r = c(r[seq_len(RECONPLOT_MAX_REASONS)], sprintf("%d more", length(r) - RECONPLOT_MAX_REASONS))
          paste(r, collapse = " · ")
        }), by = chrom]
}

# Slim mode embeds the per-chromosome figures of at most this many chromosomes (per combination),
# best reason first: chromothripsis (lower tier first), Padfoot, tier-1 then tier-2 findings
RECONPLOT_MAX_CHROMS = 12L

# reconplot_chromosomes() rows ranked and split at `n`: list(kept, dropped), both in rank order
reconplot_cap = function(sel, n = RECONPLOT_MAX_CHROMS, chromosomes = NULL) {
  if (is.null(sel) || nrow(sel) == 0) return(list(kept = sel, dropped = sel))
  cr = if (!is.null(chromosomes)) match(sel$chrom, chromosomes) else seq_len(nrow(sel))
  o = order(sel$priority, is.na(sel$ss_tier), sel$ss_tier, is.na(cr), cr)
  r = sel[o]
  k = seq_len(min(n, nrow(r)))
  list(kept = r[k], dropped = r[-k])
}

RECONPLOT_COMBOS = c("severus_ascat", "severus_wakhan", "savana")
RECONPLOT_COMBO_LABELS = c(severus_ascat = "Severus + ASCAT", severus_wakhan = "Severus + Wakhan",
                           savana = "SAVANA")

register_section(list(
  id     = "reconplot",
  title  = "ReConPlot",

  locate = function(sample_dir, sample_id) {
    # Recursive file listing sees the staged symlinked files (not symlinked directories)
    pngs = list.files(sample_dir, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
    m = regmatches(pngs, regexec("(^|/)reconplot/([^/]+)/(genome_wide|per_chromosome|focus)/([^/]+)\\.png$",
                                 pngs))
    ok = lengths(m) > 0
    if (!any(ok)) return(list(files = data.table(combo = character(), kind = character(),
                                                 file = character(), path = character()),
                              sample_id = sample_id))
    parts = do.call(rbind, m[ok])
    list(files = data.table(combo = parts[, 3], kind = parts[, 4], file = parts[, 5],
                            path = pngs[ok]),
         sample_id = sample_id)
  },

  parse  = function(inputs, section_data, ctx) {
    f = inputs$files
    if (is.null(f) || nrow(f) == 0) return(NULL)
    sid = inputs$sample_id

    # Several samples' figures in one tree: keep this sample's, where it can be told apart
    own = function(d, suffix) {
      hit = d[file == paste0(sid, suffix)]
      if (nrow(hit) > 0) hit else d
    }

    combos = unique(f$combo)
    combos = c(intersect(RECONPLOT_COMBOS, combos), sort(setdiff(combos, RECONPLOT_COMBOS)))
    chrom_order = if (!is.null(ctx$chromosomes)) as.character(unlist(ctx$chromosomes)) else character(0)

    out = list()
    for (cb in combos) {
      g = own(f[combo == cb & kind == "genome_wide" & grepl("_genome_wide$", file)],
              "_genome_wide")
      genome = if (nrow(g) > 0) g$path[order(nchar(g$path))][1] else NULL

      pc = f[combo == cb & kind == "per_chromosome" & grepl("_chr[^_]+$", file)]
      if (nrow(pc) > 0) {
        pc[, chrom := sub("^.*_(chr[^_]+)$", "\\1", file)]
        mine = pc[startsWith(file, paste0(sid, "_chr"))]
        if (nrow(mine) > 0) pc = mine
        pc = pc[!duplicated(chrom)]
        rank = match(pc$chrom, chrom_order)
        # Off-list contigs after the known ones, in natural order (chr2 before chr10)
        nat = suppressWarnings(as.numeric(sub("^chr", "", pc$chrom)))
        pc = pc[order(is.na(rank), rank, is.na(nat), nat, chrom)]
        chroms = setNames(pc$path, pc$chrom)
      } else {
        chroms = setNames(character(0), character(0))
      }

      fc = f[combo == cb & kind == "focus"]
      if (nrow(fc) > 0) {
        mine = fc[startsWith(file, paste0(sid, "_"))]
        if (nrow(mine) > 0) fc = mine
        fc = fc[order(file)]
        lab = ifelse(startsWith(fc$file, paste0(sid, "_")), substring(fc$file, nchar(sid) + 2L), fc$file)
        focus = setNames(fc$path, gsub("_", " ", lab, fixed = TRUE))
      } else {
        focus = setNames(character(0), character(0))
      }

      if (!is.null(genome) || length(chroms) > 0 || length(focus) > 0)
        out[[cb]] = list(genome = genome, chroms = chroms, focus = focus)
    }
    if (length(out) == 0) return(NULL)
    list(combos = out)
  }
))
