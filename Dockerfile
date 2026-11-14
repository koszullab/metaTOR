# Built by the release workflow from the released wheel (dist/*.whl), so the
# image contains exactly what was published to PyPI.
FROM ghcr.io/mamba-org/micromamba:2.10.0

ARG PY_VERSION=3.11

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Non-Python tools: aligners, prodigal, hmmer...
COPY --chown=$MAMBA_USER:$MAMBA_USER metator.yaml /tmp/metator.yaml
RUN micromamba install -y -n base -f /tmp/metator.yaml python=$PY_VERSION \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER dist/*.whl /tmp/dist/
RUN micromamba run -n base pip install --no-cache-dir /tmp/dist/*.whl \
    && rm -rf /tmp/dist /tmp/metator.yaml

WORKDIR /home/mambauser
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
