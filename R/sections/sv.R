# Structural variants section (reference implementation of the section-module contract), keyed by caller
#
# Callers: Severus (variants/severus/.../severus_somatic.vcf.gz + vep/SVs/<id>_SV_VEP.vcf.gz) and
# SAVANA (variants/savana/<id>.[pon_filtered.]classified.somatic.vcf + vep/SVs/<id>_VEP_SAVANA.vcf.gz),
# each optionally annotated by Padfoot (padfoot/<combo>/annotated_svs.tsv).
#
# parse() returns
#   table          one row per rearrangement, both callers rbind-ed. Columns: SV_TABLE_COLS
#                  (id, svclass, svtype, chrom_a, pos_a, gene_a, chrom_b, pos_b, gene_b, sv_len,
#                  vaf, consequence, impact, support, svtype_raw, filter, wl_rescue) + caller,
#                  and PADFOOT_COLS
#                  (padfoot_impact, padfoot_score, padfoot_role, padfoot_genes, padfoot_type)
#                  only when at least one Padfoot file was found. `id` is unique across callers
#                  (a caller whose ids collide with an earlier one's is prefixed "<caller>:").
#   circos         genome-wide circos tracks (severus_circos_tracks()) from ONE caller, without
#                  its whitelist-rescued rows (the static plot shows the default view)
#   circos_caller  which one: Severus when it has records, else SAVANA, else NA
#   circos_rescued rows of the circos caller left out of the circos as whitelist-rescued
#   n_by_caller    named integer, rows per caller found (0 for a caller whose VCF was empty)
#   n_rescued_by_caller  named integer, whitelist-rescued rows among them (sv_is_rescued())
#   callers        per caller: list(vcf, vep_vcf, gene_tsv, annotation_path, n_records, n_rows,
#                  n_annotated, n_with_genes, n_typed, pon_filtered, n_rescued,
#                  n_rescued_fail, n_rescued_1read, n_failed)
#   padfoot        per caller with a Padfoot file: list(path, combo, key, n_joined, n_rows, n_padfoot)
#   has_padfoot    TRUE when any Padfoot file was found (then PADFOOT_COLS exist in a non-empty
#                  `table`, NA on the rows of a caller without one)
#   annotation_path, annotation_found   the first caller's annotation source (back-compat)
#
# Counting rearrangements over `table` counts each one once PER CALLER: for a header figure, take
# the circos_caller's rows (table[caller == circos_caller]).

# Preferred caller order: the circos and the header figures follow the first one with records
SV_CALLERS = c("severus", "savana")

register_section(list(
  id = "sv",
  title = "Structural variants",

  locate = function(sample_dir, sample_id) {
    # One recursive walk for every SV input (list.files() sees staged symlinked files)
    all = list.files(sample_dir, recursive = TRUE, full.names = TRUE)
    base = basename(all)

    # Several hits: the one named for this sample, then the shortest path
    pick = function(hits) {
      if (length(hits) == 0) return(NULL)
      own = startsWith(basename(hits), sample_id)
      if (any(own)) hits = hits[own]
      hits[order(nchar(hits), hits)][1]
    }
    find1 = function(pattern) pick(all[grepl(pattern, base)])

    severus_vcf      = find1("^severus_somatic\\.vcf\\.gz$")
    severus_gene_tsv = find1("^SV_filtered_with_gene_annotations\\.tsv$")
    # VEP SV VCF is the usual annotation source; the gene-annotated TSV is the fallback.
    # `_SV_VEP.vcf.gz` cannot match SAVANA's `_VEP_SAVANA.vcf.gz`, nor the reverse
    severus_vep_vcf  = find1("_SV_VEP\\.vcf\\.gz$")
    savana_vep_vcf   = find1("_VEP_SAVANA\\.vcf\\.gz$")

    # SAVANA's somatic calls: the PoN-filtered file where the pipeline wrote one (PacBio
    # tumour-only; it is what VEP/Padfoot/ReConPlot were fed), else the classified somatic VCF.
    # Not classified.vcf (somatic + germline), classified.germline.vcf, sv_breakpoints.vcf, or
    # a lenient/strict variant
    savana_pon = find1("\\.pon_filtered\\.classified\\.somatic\\.vcf(\\.gz)?$")
    savana_vcf = if (!is.null(savana_pon)) savana_pon else
      pick(all[grepl("\\.classified\\.somatic\\.vcf(\\.gz)?$", base) &
               !grepl("pon_filtered|lenient|strict", base)])

    # Padfoot: padfoot/<combo>/annotated_svs.tsv; Severus is annotated with Wakhan CN
    # (severus_wakhan), SAVANA with its own (savana)
    pf = all[base == "annotated_svs.tsv"]
    pf_combo = basename(dirname(pf))
    padfoot_for = function(prefer, pattern) {
      i = which(pf_combo == prefer)
      if (length(i) == 0) i = which(grepl(pattern, pf_combo))
      if (length(i) == 0) return(NULL)
      list(path = pf[i[1]], combo = pf_combo[i[1]])
    }

    list(callers = list(
      severus = list(vcf = severus_vcf, gene_tsv = severus_gene_tsv, vep_vcf = severus_vep_vcf,
                     padfoot = padfoot_for("severus_wakhan", "^severus")),
      savana  = list(vcf = savana_vcf, vep_vcf = savana_vep_vcf,
                     pon_filtered = !is.null(savana_pon),
                     padfoot = padfoot_for("savana", "^savana"))
    ))
  },

  parse = function(inputs, section_data) {
    tabs = list(); records_by = list(); prov = list(); pad = list()
    has_padfoot = any(vapply(inputs$callers, function(ci) !is.null(ci$padfoot), logical(1)))

    for (nm in names(inputs$callers)) {
      ci = inputs$callers[[nm]]
      found = !is.null(ci$vcf) || !is.null(ci$gene_tsv)
      if (!found) next

      # One parse feeds the table and the circos tracks
      records = parse_sv_vcf(ci$vcf, nm)
      n_records = attr(records, "n_records")
      if (is.null(n_records)) n_records = 0L
      annotation_path = NULL
      n_annot = 0L

      if (nrow(records) > 0) {
        # A missing VEP SV VCF only empties the symbol columns; panel matching is on coordinates
        t = build_sv_table_from_vep(records, ci$vep_vcf)
        n_annot = attr(t, "n_annotated")
        if (!is.null(ci$vep_vcf)) annotation_path = ci$vep_vcf
      } else if (!is.null(ci$gene_tsv)) {
        # No caller VCF: fall back to the gene-annotated TSV (Severus only)
        t = build_sv_table(parse_severus_gene_tsv(ci$gene_tsv))
        if (nrow(t) > 0) annotation_path = ci$gene_tsv
      } else {
        t = data.table()
      }

      # Padfoot joins on the caller's own ids, so before any cross-caller prefixing
      if (has_padfoot && nrow(t) > 0) {
        pf = if (!is.null(ci$padfoot)) parse_padfoot(ci$padfoot$path) else NULL
        j  = join_padfoot(t, pf)
        t  = j$table
        if (!is.null(ci$padfoot))
          pad[[nm]] = list(path = ci$padfoot$path, combo = ci$padfoot$combo, key = j$key,
                           n_joined = j$n_joined, n_rows = j$n_rows,
                           n_padfoot = if (is.null(pf)) 0L else nrow(pf))
      }

      if (nrow(t) > 0) {
        t[, caller := nm]
        # Ids unique within the caller (a VCF may reuse "." or an id)
        bad = is.na(t$id) | t$id %in% c("", ".")
        if (any(bad)) t[bad, id := paste0(nm, "_row", which(bad))]
        if (anyDuplicated(t$id)) t[, id := make.unique(id, sep = "_dup")]
        tabs[[nm]] = t
      }
      records_by[[nm]] = records
      prov[[nm]] = list(
        vcf = ci$vcf, vep_vcf = ci$vep_vcf, gene_tsv = ci$gene_tsv,
        annotation_path = annotation_path,
        n_records = as.integer(n_records), n_rows = nrow(t),
        n_annotated = as.integer(n_annot),
        n_with_genes = if (nrow(t) > 0) sum(!is.na(t$gene_a) | !is.na(t$gene_b)) else 0L,
        n_typed = if (nrow(t) > 0 && "svtype_raw" %in% names(t))
                    sum(t$svtype_raw %in% BND_SVTYPES & !t$svtype %in% BND_SVTYPES) else 0L,
        pon_filtered = isTRUE(ci$pon_filtered),
        n_rescued = sum(sv_is_rescued(t)),
        n_rescued_fail = if (nrow(t) > 0 && "filter" %in% names(t))
                           sum(sv_is_rescued(t) & !is.na(t$filter) & !t$filter %in% SV_PASS_FILTERS) else 0L,
        n_rescued_1read = if (nrow(t) > 0 && "support" %in% names(t))
                            sum(sv_is_rescued(t) & t$support %in% 1L) else 0L,
        n_failed = { x = attr(records, "n_failed"); if (is.null(x)) 0L else as.integer(x) }
      )
    }

    # Ids unique across callers too: the table's row selection and the breakend circos key on
    # them. A colliding caller is prefixed as a whole, so its ids stay consistently named
    seen = character(0)
    for (nm in names(tabs)) {
      if (any(tabs[[nm]]$id %in% seen)) tabs[[nm]][, id := paste0(nm, ":", id)]
      seen = c(seen, tabs[[nm]]$id)
    }

    tbl = if (length(tabs) > 0) rbindlist(tabs, use.names = TRUE, fill = TRUE) else data.table()

    # The genome-wide circos draws one caller, never both overlaid
    with_records = names(records_by)[vapply(records_by, nrow, integer(1)) > 0]
    circos_caller = intersect(c(SV_CALLERS, setdiff(with_records, SV_CALLERS)), with_records)
    circos_caller = if (length(circos_caller)) circos_caller[1] else NA_character_
    circ_rec = if (is.na(circos_caller)) NULL else records_by[[circos_caller]]
    circ_resc = sv_is_rescued(circ_rec)
    circ = severus_circos_tracks(if (is.null(circ_rec)) NULL else circ_rec[!circ_resc])

    n_by_caller = vapply(prov, function(p) as.integer(p$n_rows), integer(1))
    n_rescued_by_caller = vapply(prov, function(p) as.integer(p$n_rescued), integer(1))

    first_annot = Filter(Negate(is.null), lapply(prov, `[[`, "annotation_path"))
    annotation_path = if (length(first_annot)) first_annot[[1]] else NULL

    list(
      table            = tbl,
      circos           = circ,
      circos_caller    = circos_caller,
      circos_rescued   = sum(circ_resc),
      n_by_caller      = n_by_caller,
      n_rescued_by_caller = n_rescued_by_caller,
      callers          = prov,
      padfoot          = pad,
      has_padfoot      = has_padfoot,
      annotation_path  = annotation_path,
      annotation_found = !is.null(annotation_path)
    )
  }
))
