# lrsomatic_report

Standalone reporting tool for the [LRSomatic](https://github.com/nf-core/lrsomatic) Nextflow pipeline. Generates a self-contained HTML report per sample with:

- **Summary header**: purity, ploidy, coverage, N50, variant counts
- **Circos plot**: somatic SNVs (6-class SBS colours), non-BND SVs, ASCAT copy number, translocation links
- **Interactive variant table**: VEP-annotated somatic small variants filtered to a gene panel of interest, with a per-caller VAF column (ClairS-TO / ClairS)
- **Interactive SV table**: Severus structural variants annotated with gene overlaps, filtered to the gene panel
- **QC details**: mosdepth coverage, samtools flagstat, cramino read stats

## Quick start

```bash
Rscript bin/render_report.R \
  --sample-dir  /path/to/DLBCL3_pooled \
  --sample-id   DLBCL3_pooled \
  --sex         male \
  --mode        matched \
  --somatic-vcf /path/to/DLBCL3_pooled/variants/clairs/somatic.vcf.gz \
  --reference   auto            # auto-detects t2t vs hg38 from VCF headers
```

The output file `DLBCL3_pooled_report.html` will be written to the current directory.

## All options

```
--sample-dir   Path to the sample output directory (required)
--sample-id    Sample identifier (default: directory name)
--reference    t2t | hg38 | auto  (default: auto)
--sex          male | female | XY | XX  (required)
--mode         matched | tumour-only  (required)
--somatic-vcf  Path to the somatic small-variant caller VCF used for VAF (required)
--gene-panel   Builtin panel name (e.g. lymphoid) or path to a custom TSV  (default: lymphoid)
--output       Output HTML path  (default: <sample-id>_report.html in current dir)
--title        Report title
```

## Gene panels

Built-in panels live in `assets/gene_lists/`. Each is a TSV with a `gene` column (HGNC symbols).

| Panel | Description |
|---|---|
| `lymphoid` | ~70 recurrently mutated genes in B-cell lymphomas (DLBCL, FL, CLL, MCL, BL, MALT) |

To use a custom panel:

```bash
--gene-panel /path/to/my_genes.tsv   # must have a 'gene' column or be a single-column file
```

## Expected input layout

The `--sample-dir` must be the root of a single-sample LRSomatic output. Files are discovered
**recursively** by their distinctive filename suffix, so they can be nested in any directory
structure underneath it — for example:

```
DLBCL3_pooled/
├── *_SOMATIC_VEP.vcf.gz                                 VEP-annotated somatic small variants
├── severus_somatic.vcf.gz                               Severus SV calls
├── SV_filtered_with_gene_annotations.tsv                Severus SV gene annotations
├── *.segments_raw.txt, *.purityploidy.txt               ASCAT
├── *.mosdepth.summary.txt, *.mosdepth.global.dist.txt    mosdepth (tumor)
├── *_cramino.txt, *.flagstat, *.stats                    cramino / samtools (tumor)
└── normal/                                               same QC file set, normal side (matched mode only)
```

The somatic small-variant caller VCF used for VAF (ClairS-TO / ClairS) has an ambiguous, generic
filename and can't be discovered reliably, so it's passed explicitly via `--somatic-vcf`; which
column it populates is determined by `--mode`.

Missing files are handled gracefully: the corresponding report section shows a "not available" notice.

## Supported references

| `--reference` | Cytobands source | chr1 length |
|---|---|---|
| `t2t` | CHM13v2.0 | 248,387,328 bp |
| `hg38` | GRCh38 (UCSC) | 248,956,422 bp |

Auto-detection reads `##contig` lines from the VEP somatic VCF.

## R package requirements

Install in your R environment if missing:

```r
install.packages(c("data.table", "dplyr", "tidyr", "DT", "htmltools",
                   "optparse", "quarto", "yaml", "ggplot2", "svglite"))
BiocManager::install(c("circlize", "ComplexHeatmap", "GenomicRanges"))
# paletteer, prismatic are optional (not required by this version)
```

Tested with R 4.4.1 and Quarto 1.5.57.

## Repository structure

```
lrsomatic_report/
├── bin/render_report.R          CLI entrypoint
├── R/
│   ├── utils.R                  Shared helpers (gene panel, Extra-field parser)
│   ├── references.R             Cytoband + chrom-length loading, reference auto-detection
│   ├── locate_outputs.R         Discover per-tool output files in a sample directory
│   ├── parse_smallvariants.R    VEP text + raw caller VCF parsers; build variant table
│   ├── parse_severus.R          Severus VCF + gene TSV parsers; build SV table
│   ├── parse_ascat.R            ASCAT segments + purity/ploidy parsers
│   ├── parse_qc.R               Mosdepth, cramino, flagstat parsers
│   └── circos.R                 draw_circos() — generates the circos SVG
├── templates/per_sample.qmd    Quarto template (HTML report)
├── assets/
│   ├── references/{t2t,hg38}/  Cytobands + chrom lengths (bundled, no network needed)
│   └── gene_lists/             lymphoid.tsv + README
└── tests/                      Unit tests (testthat)
```

## Roadmap

- **v2**: Cohort report (oncoprint, recurrence tables across multiple samples)
- **v2**: Nextflow module wrapping this CLI as a final pipeline step
- **v2**: Wakhan haplotype-resolved copy-number integration
