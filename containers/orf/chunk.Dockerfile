# ==============================================================================
# Dockerfile.chunk - Minimal container for 'orf chunk' subcommand
# ==============================================================================
# This is a standalone Rust binary with no external dependencies.
# Result: ~20-30MB container

# ---------- Build Stage ----------
FROM rust:1.93-bullseye AS builder

WORKDIR /build

COPY modules/orf/Cargo.toml modules/orf/Cargo.lock ./
COPY modules/orf/src ./src

RUN cargo build --release && \
    strip target/release/orf

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim

# Install minimal runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy the binary
COPY --from=builder /build/target/release/orf /usr/local/bin/orf

# Set up non-root user
RUN useradd -m -u 1000 orfuser && \
    chmod +x /usr/local/bin/orf

USER orfuser
WORKDIR /data

# Test that it works
RUN orf chunk --help
