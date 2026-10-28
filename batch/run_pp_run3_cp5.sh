#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 JOB_ID EVENTS NCH_MIN CAMPAIGN_SEED EOS_DIR" >&2
  exit 2
fi

job_id="$1"
events="$2"
nch_min="$3"
campaign_seed="$4"
eos_dir="$5"

output_basename="pp_run3_cp5_13p6_nchGE${nch_min}_job${job_id}.root"
eos_output="root://eoscms.cern.ch/${eos_dir}/${output_basename}"

source ./setup_lxplus.sh

# Ensure EOS output directory exists
xrdfs root://eoscms.cern.ch mkdir -p "${eos_dir}" || true

./pp_run3_cp5_13p6 \
  --job-id        "${job_id}" \
  --events        "${events}" \
  --nch-min       "${nch_min}" \
  --campaign-seed "${campaign_seed}" \
  --output        "${eos_output}"

echo "Done: ${eos_output}"
