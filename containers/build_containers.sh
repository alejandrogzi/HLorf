#!/usr/bin/env bash
# ==============================================================================
# build_containers.sh - Build all pipeline containers
# ==============================================================================
# 
# USAGE:
#   Place this script in: <project_root>/containers/
#   Run from project root: ./containers/build_containers.sh
#   OR run from script dir: cd containers && ./build_containers.sh
#
# REQUIREMENTS:
#   - Dockerfiles must be in: containers/orf/Dockerfile.* and containers/predict/Dockerfile.*
#   - Rust source must be in: modules/orf/
#   - Python sources must be in: modules/orf/tai/ and modules/predict/
#
# project_root/
# ├── containers/
# │   ├── build_containers.sh  ← This script
# │   ├── orf/
# │   │   ├── Dockerfile.chunk
# │   │   ├── Dockerfile.tai
# │   │   ├── Dockerfile.samba
# │   │   └── Dockerfile.net
# │   │   └── Dockerfile.blast
# │   └── predict/
# │       └── Dockerfile.predict
# └── modules/
#     ├── orf/
#     └── predict/
#
# ==============================================================================

set -euo pipefail

# Detect project root (should contain modules/ directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if we're in containers/ and go up to project root
if [[ "${SCRIPT_DIR}" == */containers ]]; then
    PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
elif [[ -d "${SCRIPT_DIR}/modules" ]]; then
    # Already in project root
    PROJECT_ROOT="${SCRIPT_DIR}"
else
    echo "ERROR: Cannot determine project root."
    echo "Expected directory structure:"
    echo "  project_root/"
    echo "  ├── modules/"
    echo "  │   ├── orf/"
    echo "  │   ├── nets/"
    echo "  │   └── predict/"
    echo "  └── containers/"
    echo "      ├── orf/"
    echo "      └── predict/"
    echo ""
    echo "Current directory: ${SCRIPT_DIR}"
    exit 1
fi

# Verify required directories exist
if [[ ! -d "${PROJECT_ROOT}/modules/orf" ]]; then
    echo "ERROR: modules/orf/ not found in ${PROJECT_ROOT}"
    exit 1
fi

if [[ ! -d "${PROJECT_ROOT}/modules/predict" ]]; then
    echo "ERROR: modules/predict/ not found in ${PROJECT_ROOT}"
    exit 1
fi

# Configuration
REGISTRY="${DOCKER_REGISTRY:-localhost}"
TAG="${DOCKER_TAG:-latest}"

echo "=============================================="
echo "Building Pipeline Containers"
echo "=============================================="
echo "Project Root:    ${PROJECT_ROOT}"
echo "Registry:        ${REGISTRY}"
echo "Tag:             ${TAG}"
echo ""

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# Track build results
declare -a SUCCESSFUL_BUILDS=()
declare -a FAILED_BUILDS=()

# Function to build a container
build_container() {
    local category=$1  # e.g., "orf" or "predict"
    local name=$2      # e.g., "chunk", "tai", "predict"
    local dockerfile=$3
    local image="${REGISTRY}/orf-${name}:${TAG}"
    local dockerfile_path="${PROJECT_ROOT}/containers/${category}/${dockerfile}"
    
    echo ""
    echo -e "${BLUE}[BUILD]${NC} Building ${category}-${name}..."
    echo "  Dockerfile: containers/${category}/${dockerfile}"
    echo "  Context:    ${PROJECT_ROOT}"
    echo "  Image:      ${image}"
    
    if [[ ! -f "${dockerfile_path}" ]]; then
        echo -e "${RED}[ERROR]${NC} Dockerfile not found: ${dockerfile_path}"
        FAILED_BUILDS+=("${category}-${name}")
        return 1
    fi
    
    if docker build \
        -f "${dockerfile_path}" \
        -t "${image}" \
        "${PROJECT_ROOT}"; then
        echo -e "${GREEN}[SUCCESS]${NC} Built ${image}"
        SUCCESSFUL_BUILDS+=("${image}")
        return 0
    else
        echo -e "${RED}[FAILED]${NC} Failed to build ${image}"
        FAILED_BUILDS+=("${category}-${name}")
        return 1
    fi
}

# Build ORF containers
echo ""
echo -e "${YELLOW}=== Building ORF Containers ===${NC}"
build_container "orf" "chunk" "chunk.Dockerfile" || true
build_container "orf" "tai" "tai.Dockerfile" || true
build_container "orf" "samba" "samba.Dockerfile" || true
build_container "orf" "blast" "blast.Dockerfile" || true
build_container "orf" "net" "net.Dockerfile" || true

# Build Predict containers
echo ""
echo -e "${YELLOW}=== Building Predict Containers ===${NC}"
build_container "predict" "predict" "Dockerfile.predict" || true

echo ""
echo "=============================================="
echo "Build Summary"
echo "=============================================="

# Show successful builds
if [[ ${#SUCCESSFUL_BUILDS[@]} -gt 0 ]]; then
    echo -e "${GREEN}Successfully built (${#SUCCESSFUL_BUILDS[@]}):${NC}"
    for image in "${SUCCESSFUL_BUILDS[@]}"; do
        echo "  ✓ ${image}"
    done
    echo ""
    echo "Container images:"
    docker images | head -n1
    docker images | grep -E "orf-" | grep "${TAG}" || echo "  (No matching images found)"
fi

# Show failed builds
if [[ ${#FAILED_BUILDS[@]} -gt 0 ]]; then
    echo ""
    echo -e "${RED}Failed builds (${#FAILED_BUILDS[@]}):${NC}"
    for name in "${FAILED_BUILDS[@]}"; do
        echo "  ✗ ${name}"
    done
    echo ""
    echo -e "${RED}Build completed with errors!${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}All containers built successfully!${NC}"
echo ""
echo "Next steps:"
echo "  1. Test containers: ./containers/test_containers.sh"
echo "  2. Push to registry: docker push ${REGISTRY}/orf-chunk:${TAG}"
echo "  3. Run pipeline: nextflow run src/main.nf"

