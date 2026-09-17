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

# The real thing: render a report as a non-root uid with HOME and TMPDIR pointed at the work
# dir, which is the configuration a Nextflow task runs in. This is what catches a Quarto or
# Deno that cannot find its own tooling or cannot write where it expects -- `--help` exits
# long before any of that is touched.
echo "== render a report as a non-root uid =="
rm -rf render && mkdir -p render/sample_dir
cp tests/fixtures/test_SOMATIC_VEP.vcf.gz render/sample_dir/
docker run --rm --user "$(id -u):$(id -g)" \
    -v "$PWD/render:/work" -w /work \
    -e HOME=/work -e TMPDIR=/work \
    "$image" render_report.R \
        --sample-dir sample_dir \
        --sample-id smoke \
        --sex male \
        --reference auto \
        --output smoke_report.html
test -s render/smoke_report.html
# The fixture carries a TP53 variant; its absence would mean an empty report rendered cleanly
grep -q TP53 render/smoke_report.html

# --gene-lists-dir must resolve a panel directory bind-mounted from outside the image, which
# is exactly how the pipeline supplies its own panels. A panel named here is a *builtin*, so
# it has to reach DEFAULT_PANELS in the rendered HTML.
echo "== --gene-lists-dir reads a mounted directory =="
rm -rf panels && mkdir -p panels
printf 'gene\nTP53\nMYC\n' > panels/smokepanel.tsv
docker run --rm --user "$(id -u):$(id -g)" \
    -v "$PWD/render:/work" -v "$PWD/panels:/panels" -w /work \
    -e HOME=/work -e TMPDIR=/work \
    "$image" render_report.R \
        --sample-dir sample_dir \
        --sample-id smoke2 \
        --sex male \
        --reference auto \
        --gene-lists-dir /panels \
        --gene-panel smokepanel \
        --output smoke2_report.html
grep -q 'const DEFAULT_PANELS' render/smoke2_report.html
grep -q 'smokepanel' render/smoke2_report.html
# That the bundled set is replaced rather than merged is asserted in tests/testthat, where
# it can be checked precisely instead of by grepping rendered HTML.

echo "all smoke tests passed"
