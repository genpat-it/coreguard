#!/bin/bash
# Generate demo_config.yaml from template by expanding environment variables
#
# Usage:
#   export SNPS_STUDY=/path/to/snps-study
#   export COREGUARD=/path/to/coreguard
#   ./scripts/make_demo_config.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Default paths (adjust for your environment)
: "${COREGUARD:=$PROJECT_DIR}"
: "${SNPS_STUDY:=$HOME/snps-study}"

TEMPLATE="$PROJECT_DIR/testdata/demo_config.template.yaml"
OUTPUT="$PROJECT_DIR/testdata/demo_config.yaml"

if [ ! -f "$TEMPLATE" ]; then
    echo "Error: Template not found: $TEMPLATE"
    exit 1
fi

echo "Generating demo_config.yaml..."
echo "  COREGUARD=$COREGUARD"
echo "  SNPS_STUDY=$SNPS_STUDY"

# Expand variables in template
envsubst '${COREGUARD} ${SNPS_STUDY}' < "$TEMPLATE" > "$OUTPUT"

echo "Created: $OUTPUT"
