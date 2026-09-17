# lrsomatic_report packaged as a container, so consumers (the LRSomatic Nextflow pipeline)
# can pull the tool instead of vendoring a copy of this source tree.
#
# The install layout below is the one recipe/build.sh produces, so once `lrsomatic-report`
# reaches bioconda the biocontainer is a drop-in replacement and callers change nothing but
# the image reference.
#
# Built by .github/workflows/container.yml on a v* tag and published as
#   ghcr.io/ljwharbers/lrsomatic-report:<version>            (OCI)
#   oras://ghcr.io/ljwharbers/lrsomatic-report-sif:<version> (Apptainer SIF)
FROM condaforge/miniforge3:26.7.2-0

ARG VERSION=dev
LABEL org.opencontainers.image.title="lrsomatic_report" \
      org.opencontainers.image.description="R/Quarto reporting tool for the LRSomatic long-read somatic variant calling pipeline" \
      org.opencontainers.image.source="https://github.com/ljwharbers/lrsomatic_report" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.licenses="MIT"

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

# procps-ng: Nextflow's task monitor calls `ps` inside the container. It is not a dependency
# of the tool, which is why it is here and not in recipe/meta.yaml.
COPY container/environment.yml /tmp/environment.yml
RUN mamba env update -n base -f /tmp/environment.yml \
 && mamba install -y -n base conda-forge::procps-ng \
 && mamba clean -afy \
 && rm -f /tmp/environment.yml

# quarto's conda package exports QUARTO_SHARE_PATH from etc/conda/activate.d/*.sh. Neither
# `apptainer exec` nor a Nextflow task sources activation hooks, so bake it in rather than
# making every caller replicate the activation dance.
ENV PATH=/opt/conda/bin:$PATH \
    QUARTO_SHARE_PATH=/opt/conda/share/quarto \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# Same layout as recipe/build.sh. render_report.R resolves its own bin/ symlink with
# normalizePath() before taking dirname()/.., so repo_dir lands on share/lrsomatic_report.
COPY bin/       /opt/conda/share/lrsomatic_report/bin/
COPY R/         /opt/conda/share/lrsomatic_report/R/
COPY templates/ /opt/conda/share/lrsomatic_report/templates/
COPY assets/    /opt/conda/share/lrsomatic_report/assets/
COPY VERSION LICENSE /opt/conda/share/lrsomatic_report/

RUN chmod +x /opt/conda/share/lrsomatic_report/bin/render_report.R \
 && ln -s ../share/lrsomatic_report/bin/render_report.R /opt/conda/bin/render_report.R \
 && chmod -R a+rX /opt/conda/share/lrsomatic_report

# Build-time assertions. Each one stands for a way this image has a plausible path to being
# silently broken: a QUARTO_SHARE_PATH that points nowhere, a missing R dependency that only
# shows up mid-render, and a bin/ symlink whose repo_dir does not resolve.
RUN test -d "$QUARTO_SHARE_PATH" \
 && quarto --version \
 && Rscript -e 'for (p in c("optparse","quarto","yaml","data.table","dplyr","DT","htmltools","ggplot2","svglite","circlize","knitr","R.utils","base64enc")) library(p, character.only = TRUE)' \
 && test "$(render_report.R --version)" = "$(cat /opt/conda/share/lrsomatic_report/VERSION)"

CMD ["render_report.R", "--help"]
