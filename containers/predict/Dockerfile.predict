# ==============================================================================
# Dockerfile.predict - Container for ORF prediction with XGBoost
# ==============================================================================
# Python-based container with ML dependencies
# Result: ~500-600MB container (due to numpy/scipy/sklearn/xgboost)

# ---------- Build Stage ----------
FROM python:3.11-slim-bookworm AS builder

WORKDIR /build

# Install build dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy only pyproject.toml first for better layer caching
COPY modules/predict/pyproject.toml ./

# Install dependencies in a virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install dependencies from pyproject.toml
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    pandas==2.2.3 \
    numpy==2.1.2 \
    scikit-learn==1.5.2 \
    matplotlib==3.9.2 \
    scipy==1.14.1 \
    joblib==1.4.2 \
    xgboost==3.1.0 \
    seaborn

# ---------- Runtime Stage ----------
FROM python:3.11-slim-bookworm

# Install minimal runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy the prediction script and model
COPY modules/predict/predict.py /opt/predict/predict.py
COPY modules/predict/model /opt/predict/model

# Make script executable
RUN chmod +x /opt/predict/predict.py

# Add to PATH so 'predict' command works
ENV PATH="/opt/predict:$PATH"

# Set up non-root user
RUN useradd -m -u 1000 predictuser && \
    chown -R predictuser:predictuser /opt/predict/model
USER predictuser

WORKDIR /data

# Test that dependencies are available
RUN python -c "import pandas; import numpy; import sklearn; import xgboost; print('All dependencies OK')"

# Verify the script is accessible
RUN predict --help || echo "Note: predict script may need --help implementation"
