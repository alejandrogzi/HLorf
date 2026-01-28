# ---------- Rust build stage ----------
FROM rust:1.93-bullseye AS builder
WORKDIR /build

COPY modules/orf/Cargo.toml modules/orf/Cargo.lock ./
COPY modules/orf/src ./src
RUN cargo build --release


# ---------- Runtime stage ----------
FROM mambaorg/micromamba:1.5.8

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
ENV MAMBA_ROOT_PREFIX=/opt/conda
WORKDIR /opt

# 1) RNAsamba env: pinned pip stack (mirrors upstream Dockerfile)
RUN micromamba create -y -n rnasambaenv -c conda-forge python=3.6 pip \
 && micromamba clean -a -y

RUN micromamba run -n rnasambaenv python -m pip install --upgrade pip "setuptools<60" wheel \
 && micromamba run -n rnasambaenv python -m pip install --no-cache-dir \
      "biopython==1.74" \
      "keras==2.2.5" \
      "numpy==1.16.5" \
      "h5py==2.10.0" \
      "tensorflow==1.14.0" \
      "rnasamba==0.2.5"

# 2) Diamond env
RUN micromamba create -y -n diamondenv -c conda-forge -c bioconda \
      python=3.9 diamond \
 && micromamba clean -a -y

# 3) translationai env (Python 3.9 + deps)
RUN micromamba create -y -n taienv -c conda-forge \
      python=3.9 pip \
      gcc_linux-64 gxx_linux-64 \
      make \
      pkg-config \
 && micromamba clean -a -y

# 4) Copy Rust binary
COPY --from=builder /build/target/release/orf /usr/local/bin/orf


USER root 

# 5) Copy translationai source (layout matches where=["./tai"])
COPY modules/orf/pyproject.toml /opt/translationai/pyproject.toml
COPY modules/orf/tai /opt/translationai/tai

# Make source writable for the default non-root user (uid 1000)
RUN chown -R 1000:1000 /opt/translationai

USER 1000
ENV HOME=/tmp
ENV XDG_CACHE_HOME=/tmp/.cache
RUN mkdir -p /tmp/.cache

# 6) Install translationai deps + translationai itself INTO taienv
RUN micromamba run -n taienv python -m pip install --upgrade pip "setuptools>=61.0" wheel \
 && micromamba run -n taienv python -m pip install \
      "numpy>=1.14.0,<2" \
      "pandas>=0.24.2,<3" \
      "h5py" \
      "scikit-learn>=0.20.0,<2" \
      "joblib>=1.0.0,<2" \
      "orfipy>=0.0.4" \
      "keras>=2.0.5,<3" \
      "tensorflow==2.10.*" \
 && micromamba run -n taienv python -m pip install /opt/translationai 

USER root

# 7) Create direct-call wrappers (no need to micromamba run manually)
#    These wrappers ensure the correct environment is used regardless of PATH.

RUN cat <<'EOF' > /usr/local/bin/rnasamba
#!/usr/bin/env bash
set -euo pipefail
exec micromamba run -n rnasambaenv rnasamba "$@"
EOF

RUN chmod 0755 /usr/local/bin/rnasamba


RUN cat <<'EOF' > /usr/local/bin/diamond
#!/usr/bin/env bash
set -euo pipefail
exec micromamba run -n diamondenv diamond "$@"
EOF

RUN chmod 0755 /usr/local/bin/diamond


RUN cat <<'EOF' > /usr/local/bin/translationai
#!/usr/bin/env bash
set -euo pipefail
exec micromamba run -n taienv translationai "$@"
EOF

RUN chmod 0755 /usr/local/bin/translationai

USER 1000

# 8) Sanity checks using the direct commands (what you wanted)
RUN orf --help >/dev/null
RUN diamond --version 
# RUN translationai --help >/dev/null
# RUN rnasamba >/dev/null

ENTRYPOINT ["orf"]
