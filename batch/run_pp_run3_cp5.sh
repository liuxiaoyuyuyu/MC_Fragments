#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 JOB_ID EVENTS NCH_MIN CAMPAIGN_SEED OUTPUT EOS_DIR" >&2
  exit 2
fi

job_id="$1"
events="$2"
nch_min="$3"
campaign_seed="$4"
output_file="$5"
eos_dir="$6"

source ./setup_lxplus.sh

./pp_run3_cp5_13p6 \
  --job-id        "${job_id}" \
  --events        "${events}" \
  --nch-min       "${nch_min}" \
  --campaign-seed "${campaign_seed}" \
  --output        "${output_file}"

# Copy output to EOS via XRootD
echo "Copying ${output_file} to ${eos_dir}/"
xrdcp -f "${output_file}" "root://eoscms.cern.ch/${eos_dir}/${output_file}"
echo "Done."
