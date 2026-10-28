#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "Usage: $0 JOB_ID EVENTS NCH_MIN CAMPAIGN_SEED SAVE_CONSTITUENTS OUTPUT EOS_DIR" >&2
  exit 2
fi

job_id="$1"
events="$2"
nch_min="$3"
campaign_seed="$4"
save_constituents="$5"
output_file="$6"
eos_dir="$7"

source ./setup_lxplus.sh

arguments=(
  --job-id        "${job_id}"
  --events        "${events}"
  --nch-min       "${nch_min}"
  --campaign-seed "${campaign_seed}"
  --output        "${output_file}"
)
if [[ "${save_constituents}" == "1" ]]; then
  arguments+=(--save-constituents)
fi

./pp_run3_cp5_13p6 "${arguments[@]}"

# Copy output to EOS via XRootD
echo "Copying ${output_file} to ${eos_dir}/"
xrdcp -f "${output_file}" "root://eoscms.cern.ch/${eos_dir}/${output_file}"
echo "Done."
