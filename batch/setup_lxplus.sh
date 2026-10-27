#!/usr/bin/env bash

# Source this file on lxplus and in every HTCondor worker job. Override
# LCG_VIEW before sourcing if your lxplus architecture uses a different view.
set -eo pipefail

LCG_VIEW="${LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt}"
if [[ ! -r "${LCG_VIEW}/setup.sh" ]]; then
  echo "LCG view not found: ${LCG_VIEW}" >&2
  echo "Set LCG_VIEW to an available ROOT+PYTHIA8+FastJet+LHAPDF view." >&2
  return 1 2>/dev/null || exit 1
fi

# shellcheck disable=SC1090
set +u
source "${LCG_VIEW}/setup.sh"
set -u

for command_name in pythia8-config fastjet-config root-config lhapdf-config; do
  if ! command -v "${command_name}" >/dev/null 2>&1; then
    echo "Missing ${command_name} after sourcing ${LCG_VIEW}/setup.sh" >&2
    return 1 2>/dev/null || exit 1
  fi
done

if command -v lhapdf >/dev/null 2>&1; then
  if ! lhapdf show NNPDF31_nnlo_as_0118 >/dev/null 2>&1; then
    echo "LHAPDF set NNPDF31_nnlo_as_0118 is not visible in this environment." >&2
    echo "Check LHAPDF_DATA_PATH or choose an LCG view containing the set." >&2
    return 1 2>/dev/null || exit 1
  fi
fi
