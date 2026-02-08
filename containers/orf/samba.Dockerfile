# ==============================================================================
# Dockerfile.samba - Container for 'orf samba' subcommand
# ==============================================================================
# Includes: orf binary + rnasamba (Python 3.6, TensorFlow 1.14)
# Result: ~1.2GB container
# Note: Uses legacy Python 3.6 stack - this is the oldest/most fragile container

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

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create necessary directories
RUN mkdir -p /opt/conda && \
    chmod -R 0777 /opt/conda

USER $MAMBA_USER

# Set environment variables
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/bin:$PATH

# Create Python 3.6 environment for rnasamba
# Pin exact versions that rnasamba expects
RUN micromamba create -y -n sambaenv -c conda-forge \
    python=3.6 \
    pip \
    && micromamba clean -a -y

# Install rnasamba with its exact dependency stack
# This follows the upstream rnasamba Docker setup
RUN micromamba run -n sambaenv pip install --no-cache-dir \
    "setuptools<60" \
    "wheel" \
    && micromamba run -n sambaenv pip install --no-cache-dir \
    "biopython==1.74" \
    "numpy==1.16.5" \
    "h5py==2.10.0" \
    "keras==2.2.5" \
    "tensorflow==1.14.0" \
    "rnasamba==0.2.5"

# Copy Rust binary
USER root
COPY --from=builder /build/target/release/orf /usr/local/bin/orf
RUN chmod +x /usr/local/bin/orf

# Set environment variables so conda environment is available
ENV PATH="/opt/conda/envs/sambaenv/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/envs/sambaenv/lib:$LD_LIBRARY_PATH"
ENV PYTHONPATH="/opt/conda/envs/sambaenv/lib/python3.6/site-packages:$PYTHONPATH"

# Create non-root user
RUN useradd -m -u 1000 orfuser
USER orfuser
WORKDIR /data

# Verify installations (tools should be directly in PATH now)
RUN orf samba --help && \
    rnasamba --help

# ENTRYPOINT ["orf"]
# CMD ["samba", "--help"]
