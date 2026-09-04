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
# ADDITIVE ONLY — no --delete. The refresh overwrites tracked code with the
# working-copy version and adds new files. It NEVER deletes anything on Drive,
# so all gitignored Drive-only content (Old Code/, Figures/, RaX figures &
# reports, baselines, *.zip, thinking/, this README) is untouched. Files you
# removed/renamed in the repo will linger as stale copies on the snapshot —
# prune those by hand if you ever want a clean mirror (deliberately not
# automated: a wrong --delete here is irreversible Drive data loss).
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
  echo ">>> APPLYING (--go). Additive only; no deletions on Drive."
fi

# No --delete by design. Excludes below only skip pointless-to-copy cruft;
# they have no destructive effect since nothing is ever removed on the receiver.
rsync -a $DRYRUN --itemize-changes \
  --exclude='.git' \
  --exclude='.DS_Store' \
  --exclude='**/.DS_Store' \
  --exclude='__pycache__/' \
  --exclude='**/__pycache__/' \
  --exclude='*.pyc' \
  --exclude='.ipynb_checkpoints/' \
  --exclude='**/.ipynb_checkpoints/' \
  --exclude='.claude/' \
  --exclude='README.SNAPSHOT.md' \
  "$SRC" "$DST"

if [[ -n "$DRYRUN" ]]; then
  echo ">>> DRY RUN complete. Re-run with --go to apply."
else
  echo ">>> Snapshot refreshed (additive; no Drive files deleted)."
fi
