FROM ghcr.io/mamba-org/micromamba:latest

ARG PY_VERSION
ARG PKG_VERSION

LABEL Name=metator Version=$PKG_VERSION

COPY --chown=$MAMBA_USER:$MAMBA_USER . ./

# Install the package
RUN micromamba install -y -n base --file metator.yaml python=$PY_VERSION && \
    micromamba run -n base pip install . --no-deps && \
    micromamba clean --all --yes && \
    rm -rf /home/mambauser/.cache && \
    rm -rf ./*

USER mambauser
WORKDIR /home/mambauser/
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
