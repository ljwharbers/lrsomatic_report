#!/usr/bin/env bash
# Smoke-test a built lrsomatic_report image.
#
#   container/smoke_test.sh <image-ref> <expected-version>
#
# Shared by the pull-request job (which builds without pushing) and the release job, so a PR
# proves the image works before a tag ever publishes one.
set -euo pipefail

image="${1:?usage: smoke_test.sh <image-ref> <expected-version>}"
version="${2:?usage: smoke_test.sh <image-ref> <expected-version>}"

echo "== version =="
docker run --rm "$image" render_report.R --version | tee v.txt
grep -qx "$version" v.txt

echo "== help lists the flags the pipeline passes =="
docker run --rm "$image" render_report.R --help > help.txt 2>&1
for flag in --sample-dir --sample-id --sex --reference --gene-panel --gene-lists-dir --output; do
    grep -q -- "$flag" help.txt || { echo "missing flag: $flag"; exit 1; }
done

# Nextflow's task monitor calls `ps` inside the container
echo "== ps =="
docker run --rm "$image" ps --version

echo "== install layout =="
docker run --rm "$image" test -f /opt/conda/share/lrsomatic_report/bin/render_report.R
docker run --rm "$image" test -d /opt/conda/share/lrsomatic_report/assets/gene_lists
docker run --rm "$image" bash -c 'test -d "$QUARTO_SHARE_PATH"'

# Render for real as a non-root uid with HOME and TMPDIR pointed at the work dir. This is the
# configuration a Nextflow task runs in, and the one that breaks when Quarto/Deno cannot write
# where they expect -- `--help` never exercises it.
echo "== quarto render as a non-root uid =="
rm -rf smoke && mkdir -p smoke
cat > smoke/smoke.qmd <<'QMD'
---
title: smoke
format: html
---

```{r}
library(ggplot2)
library(DT)
library(data.table)
1 + 1
```
QMD
docker run --rm --user "$(id -u):$(id -g)" \
    -v "$PWD/smoke:/smoke" -w /smoke \
    -e HOME=/smoke -e TMPDIR=/smoke \
    "$image" quarto render smoke.qmd --to html
test -f smoke/smoke.html

# --gene-lists-dir must resolve a panel directory bind-mounted from outside the image, which
# is exactly how the pipeline supplies its own panels.
echo "== --gene-lists-dir reads a mounted directory =="
rm -rf panels && mkdir -p panels
printf 'gene\nTP53\nMYC\n' > panels/smokepanel.tsv
docker run --rm --user "$(id -u):$(id -g)" \
    -v "$PWD/panels:/panels" -w /tmp \
    -e HOME=/tmp -e TMPDIR=/tmp \
    "$image" Rscript -e '
      source("/opt/conda/share/lrsomatic_report/R/utils.R")
      p <- load_all_gene_panels("/panels", "hg38")
      stopifnot(identical(names(p), "smokepanel"))
      stopifnot(setequal(p[["smokepanel"]]$genes, c("TP53", "MYC")))
      cat("gene-lists-dir OK\n")'

echo "all smoke tests passed"
