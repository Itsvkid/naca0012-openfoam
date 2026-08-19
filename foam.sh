#!/usr/bin/env bash
# Run an OpenFOAM command against ./case inside the container.
#
#   ./foam.sh blockMesh
#   ./foam.sh "blockMesh && checkMesh"
#
# The `cd /case` is not redundant with docker's -w. The image's entrypoint runs
# a login shell that sources a profile which cd's to the home directory, so -w
# is silently overridden and OpenFOAM reports `Case : /root` — then fails
# looking for system/controlDict there. Everything written would land in the
# container's own filesystem and vanish with --rm.
set -euo pipefail
IMAGE="${FOAM_IMAGE:-opencfd/openfoam-default:latest}"
CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/case"
exec docker run --rm -v "$CASE_DIR:/case" "$IMAGE" bash -lc "cd /case && $*"
