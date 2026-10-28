#!/usr/bin/env bash
# Source this file on lxplus and in every HTCondor worker job.
# Uses custom AFS installs of Pythia8, FastJet, and LHAPDF6.
set -eo pipefail

PYTHIA8_PREFIX=/afs/cern.ch/user/x/xiaoyul/pythia8310_install
FASTJET_PREFIX=/afs/cern.ch/user/x/xiaoyul/fastjet-install
LHAPDF_PREFIX=/afs/cern.ch/user/x/xiaoyul/LHAPDF_Lib

export PATH="${PYTHIA8_PREFIX}/bin:${FASTJET_PREFIX}/bin:${LHAPDF_PREFIX}/bin:${PATH}"
export LD_LIBRARY_PATH="${PYTHIA8_PREFIX}/lib:${FASTJET_PREFIX}/lib:${LHAPDF_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export LHAPDF_DATA_PATH="${LHAPDF_PREFIX}/share/LHAPDF"

for cmd in pythia8-config fastjet-config root-config lhapdf-config; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Missing ${cmd} — check install paths in setup_lxplus.sh" >&2
    return 1 2>/dev/null || exit 1
  fi
done

if command -v lhapdf >/dev/null 2>&1; then
  if ! lhapdf show NNPDF31_nnlo_as_0118 >/dev/null 2>&1; then
    echo "LHAPDF set NNPDF31_nnlo_as_0118 not found in LHAPDF_DATA_PATH=${LHAPDF_DATA_PATH}" >&2
    return 1 2>/dev/null || exit 1
  fi
fi
