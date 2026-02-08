# ==============================================================================
# Dockerfile.tai - Container for 'orf tai' subcommand
# ==============================================================================
# Includes: orf binary + translationai (Python 3.9, TensorFlow 2.10, models)
# Result: ~1.5GB container

# ---------- Build Stage ----------
FROM rust:1.93-slim-bullseye AS builder

WORKDIR /build

COPY modules/orf/Cargo.toml modules/orf/Cargo.lock ./
COPY modules/orf/src/ ./src/

RUN cargo build --release && \
    strip target/release/orf

# ---------- Runtime Stage ----------
FROM mambaorg/micromamba:1.5.8

USER root

# Install system dependencies including gcc for building Python packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    gcc \
    g++ \
    procps \
    make \
    && rm -rf /var/lib/apt/lists/*

# Create necessary directories with proper permissions
RUN mkdir -p /opt/translationai /opt/conda && \
    chmod -R 0777 /opt/conda /opt/translationai

USER $MAMBA_USER

# Set environment variables for micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/bin:$PATH

# Create Python 3.9 environment with scientific stack from conda (faster than pip)
RUN micromamba create -y -n taienv -c conda-forge -c bioconda \
    python=3.9 \
    pip \
    "numpy>=1.14.0,<2" \
    "pandas>=0.24.2,<3" \
    "h5py" \
    "scikit-learn>=0.20.0,<2" \
    "joblib>=1.0.0,<2" \
    && micromamba clean -a -y

# Install TensorFlow and Keras
# Using TF 2.10.x as it's the last version with native Keras compatibility for Python 3.9
RUN micromamba run -n taienv pip install --no-cache-dir \
    "tensorflow>=2.10.0,<2.11" \
    "keras>=2.10.0,<3.0"

# Copy translationai source code
USER root
COPY modules/orf/pyproject.toml /opt/translationai/pyproject.toml
COPY modules/orf/tai /opt/translationai/tai
RUN chown -R $MAMBA_USER:$MAMBA_USER /opt/translationai

# Install orfipy and translationai as MAMBA_USER (needs write access to build/)
USER $MAMBA_USER
WORKDIR /opt/translationai
RUN micromamba run -n taienv pip install --no-cache-dir \
    "orfipy>=0.0.4" \
    /opt/translationai

# Copy Rust binary
USER root
COPY --from=builder /build/target/release/orf /usr/local/bin/orf
RUN chmod +x /usr/local/bin/orf

# Instead of wrapper scripts, we'll activate the environment in the shell init
# and create symlinks to the actual binaries
RUN micromamba run -n taienv bash -c 'ln -sf $(which translationai) /usr/local/bin/translationai' && \
    micromamba run -n taienv bash -c 'ln -sf $(which orfipy) /usr/local/bin/orfipy' && \
    micromamba run -n taienv bash -c 'ln -sf $(which python) /usr/local/bin/python3.9'

# Set environment variables so conda environment is available
ENV PATH="/opt/conda/envs/taienv/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/envs/taienv/lib:$LD_LIBRARY_PATH"
ENV PYTHONPATH="/opt/conda/envs/taienv/lib/python3.9/site-packages:$PYTHONPATH"
# Set environment variables so conda environment is available
ENV PATH="/opt/conda/envs/taienv/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/envs/taienv/lib:$LD_LIBRARY_PATH"
ENV PYTHONPATH="/opt/conda/envs/taienv/lib/python3.9/site-packages:$PYTHONPATH"

# Set up environment variable for model path
ENV TRANSLATIONAI_MODEL_PATH=""

# Create non-root user
RUN useradd -m -u 1000 orfuser
USER orfuser
WORKDIR /data

# Verify installations (tools should be directly in PATH now)
RUN orf tai --help && \
    translationai --help

# ENTRYPOINT ["orf"]
# CMD ["tai", "--help"]
