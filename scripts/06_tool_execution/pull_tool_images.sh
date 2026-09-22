#!/usr/bin/env bash
# Pull the nine Docker Hub images (What the Phage image set) used by the CPU run
# scripts (run_cpu_tools.sh, run_secondary_cpu_tools.sh) into
#   ${VBENCH_ROOT}/What_the_Phage/singularity_images/
# under the exact file names those scripts expect. Existing files are skipped.
#
# Run on a machine with internet access (compute nodes often have none):
#   source config/paths.sh
#   bash scripts/06_tool_execution/pull_tool_images.sh [--with-phigaro]
#
#   --with-phigaro  also pull the Phigaro image used by scripts/04_ground_truth/run_phigaro.sh
#
# Set APPTAINER_CACHEDIR (or SINGULARITY_CACHEDIR) to a large disk if $HOME is small.
# CONTAINER_CMD overrides the runtime (default: apptainer, else singularity).
# Image list, versions and registry digests: containers/containers.tsv
set -euo pipefail

: "${VBENCH_ROOT:?is not set. Copy config/paths.example.sh to config/paths.sh, edit it, then run: source config/paths.sh}"

IMG_DIR="${VBENCH_ROOT}/What_the_Phage/singularity_images"

WITH_PHIGARO=0
for arg in "$@"; do
    case "${arg}" in
        --with-phigaro) WITH_PHIGARO=1 ;;
        -h|--help) sed -n '2,15p' "$0"; exit 0 ;;
        *) echo "ERROR: unknown option: ${arg} (see --help)" >&2; exit 2 ;;
    esac
done

if [[ -n "${CONTAINER_CMD:-}" ]]; then
    RUNTIME="${CONTAINER_CMD}"
elif command -v apptainer >/dev/null 2>&1; then
    RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    RUNTIME="singularity"
else
    echo "ERROR: neither apptainer nor singularity is in PATH" >&2
    exit 1
fi

# "<file name expected by the run scripts> <Docker Hub image>"
IMAGES=(
    "metaphinder_0.1.img multifractal/metaphinder:0.1"
    "sourmash_4.5.0--e12a57a.img nanozoo/sourmash:4.5.0--e12a57a"
    "virsorter_0.1.2.img multifractal/virsorter:0.1.2"
    "virfinder_0.2.img multifractal/virfinder:0.2"
    "virsorter-2_2.2.1--fa935f8.img papanikos/virsorter-2:2.2.1--fa935f8"
    "vibrant_0.5.img multifractal/vibrant:0.5"
    "deepvirfinder_0.1.img multifractal/deepvirfinder:0.1"
    "ppr-meta_0.3.1.img multifractal/ppr-meta:0.3.1"
    "seeker_0.1.img multifractal/seeker:0.1"
)
if (( WITH_PHIGARO )); then
    IMAGES+=("phigaro_0.5.2.img multifractal/phigaro:0.5.2")
fi

mkdir -p "${IMG_DIR}"
echo "Runtime:   ${RUNTIME}"
echo "Image dir: ${IMG_DIR}"

# Pull to a temporary name and rename on success, so an interrupted pull never
# leaves a truncated file that a later run would skip.
PARTIAL=""
cleanup() {
    if [[ -n "${PARTIAL}" && -e "${PARTIAL}" ]]; then
        rm -f "${PARTIAL}"
    fi
}
trap cleanup EXIT

for entry in "${IMAGES[@]}"; do
    read -r file ref <<< "${entry}"
    dest="${IMG_DIR}/${file}"
    if [[ -s "${dest}" ]]; then
        echo "[skip] ${file} (already present)"
        continue
    fi
    echo "[pull] docker://${ref} -> ${file}"
    PARTIAL="${dest}.partial.$$"
    rm -f "${PARTIAL}"
    if ! "${RUNTIME}" pull "${PARTIAL}" "docker://${ref}"; then
        echo "ERROR: pull failed for docker://${ref} (target ${dest})" >&2
        exit 1
    fi
    mv "${PARTIAL}" "${dest}"
    PARTIAL=""
done

# Final check: every expected image file exists and is non-empty.
missing=0
for entry in "${IMAGES[@]}"; do
    read -r file ref <<< "${entry}"
    if [[ -s "${IMG_DIR}/${file}" ]]; then
        echo "  ok       ${file}"
    else
        echo "  MISSING  ${file} (docker://${ref})" >&2
        missing=$((missing + 1))
    fi
done
if (( missing > 0 )); then
    echo "ERROR: ${missing} expected image file(s) missing in ${IMG_DIR}" >&2
    exit 1
fi
echo "All ${#IMAGES[@]} expected image files present in ${IMG_DIR}"

echo ""
echo "Not pulled by this script:"
echo "  geNomad: run bash scripts/06_tool_execution/setup_genomad.sh"
echo "           (pulls docker://quay.io/biocontainers/genomad:1.11.2--pyhdfd78af_0 to"
echo "            containers/genomad/genomad.sif and downloads the geNomad database)"
echo "  Jaeger, HVSeeker, ViraLM, TransGINmer: build from the definitions in containers/"
echo "           (see containers/README.md)"
