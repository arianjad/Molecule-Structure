#!/usr/bin/env bash
# Refresh the Google Drive git-less snapshot from this working copy.
#
# This working copy (git / GitHub `arianjad/Molecule-Structure`) is the source
# of truth. The Drive path is a read-only snapshot of known-good code with NO
# git metadata — git on the Drive CloudStorage FUSE mount corrupts the index/
# refs, which is why the two were decoupled (2026-05-18).
#
# Dry-run by default. Pass --go to actually apply.
#
#   scripts/refresh-drive-snapshot.sh        # preview (rsync -n), writes nothing
#   scripts/refresh-drive-snapshot.sh --go   # apply
#
# --delete mirrors deletions, but excluded paths are NEVER deleted on the
# receiver (no --delete-excluded), so the Drive-only heavies survive:
# Figures/ (~539 MB), thinking/ (scratch), *.zip archives, README.SNAPSHOT.md.
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/"
DST="/Users/arianjadbabaie/Library/CloudStorage/GoogleDrive-arianjad@mit.edu/Shared drives/EMA-data-server/RaX/Personal/ArianJadbabaie/Code/Molecule-Structure/"

DRYRUN="-n"
if [[ "${1:-}" == "--go" ]]; then DRYRUN=""; fi

if [[ ! -d "$DST" ]]; then
  echo "ERROR: Drive snapshot path not found (Drive not mounted / path moved):"
  echo "  $DST"
  exit 1
fi

echo "SRC (working copy):   $SRC"
echo "DST (Drive snapshot): $DST"
if [[ -n "$DRYRUN" ]]; then
  echo ">>> DRY RUN — nothing will be written. Re-run with --go to apply."
else
  echo ">>> APPLYING (--go)."
fi

rsync -a --delete $DRYRUN --itemize-changes \
  --exclude='.git' \
  --exclude='.DS_Store' \
  --exclude='**/.DS_Store' \
  --exclude='__pycache__/' \
  --exclude='**/__pycache__/' \
  --exclude='*.pyc' \
  --exclude='.ipynb_checkpoints/' \
  --exclude='**/.ipynb_checkpoints/' \
  --exclude='.claude/' \
  --exclude='Figures/' \
  --exclude='thinking/' \
  --exclude='*.zip' \
  --exclude='README.SNAPSHOT.md' \
  "$SRC" "$DST"

if [[ -n "$DRYRUN" ]]; then
  echo ">>> DRY RUN complete. Re-run with --go to apply."
else
  echo ">>> Snapshot refreshed."
fi
