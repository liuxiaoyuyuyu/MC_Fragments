#!/usr/bin/env bash
# Source this file on lxplus and in every HTCondor worker job.
# Libraries: custom AFS installs of Pythia8, FastJet, and LHAPDF6.
# PDF data:  CVMFS (available on all lxplus/condor nodes, no local install needed).
#
# LHAPDF note: lhapdf-config has a broken baked-in prefix, so we bypass it
# by exporting LHAPDF_PREFIX directly. Makefile.run3 uses this instead.
#
# set -e is turned off at the end so sourcing this file does not make the
# interactive shell exit on subsequent command failures (e.g. a failed make).

set -eo pipefail

PYTHIA8_PREFIX=/afs/cern.ch/user/x/xiaoyul/pythia8310_install
FASTJET_PREFIX=/afs/cern.ch/user/x/xiaoyul/fastjet-install
CVMFS_LHAPDF=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current

export LHAPDF_PREFIX=/afs/cern.ch/user/x/xiaoyul/LHAPDF_Lib

export PATH="${PYTHIA8_PREFIX}/bin:${FASTJET_PREFIX}/bin:${PATH}"
export LD_LIBRARY_PATH="${PYTHIA8_PREFIX}/lib:${FASTJET_PREFIX}/lib:${LHAPDF_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export LHAPDF_DATA_PATH="${CVMFS_LHAPDF}"

for cmd in pythia8-config fastjet-config root-config; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Missing ${cmd} — check install paths in setup_lxplus.sh" >&2
    set +eo pipefail
    return 1 2>/dev/null || exit 1
  fi
done

if [[ ! -f "${LHAPDF_PREFIX}/lib/libLHAPDF.so" ]]; then
  echo "libLHAPDF.so not found at ${LHAPDF_PREFIX}/lib" >&2
  set +eo pipefail
  return 1 2>/dev/null || exit 1
fi

if [[ ! -d "${CVMFS_LHAPDF}/NNPDF31_nnlo_as_0118" ]]; then
  echo "NNPDF31_nnlo_as_0118 not found at ${CVMFS_LHAPDF}" >&2
  echo "Check that CVMFS is mounted on this node." >&2
  set +eo pipefail
  return 1 2>/dev/null || exit 1
fi

set +eo pipefail
echo "Environment ready."
echo "  pythia8-config: $(which pythia8-config)"
echo "  fastjet-config: $(which fastjet-config)"
echo "  LHAPDF_PREFIX:  ${LHAPDF_PREFIX}"
echo "  LHAPDF_DATA_PATH: ${LHAPDF_DATA_PATH}"
