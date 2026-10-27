#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 JOB_ID EVENTS NCH_MIN CAMPAIGN_SEED SAVE_CONSTITUENTS OUTPUT" >&2
  exit 2
fi

job_id="$1"
events="$2"
nch_min="$3"
campaign_seed="$4"
save_constituents="$5"
output_file="$6"

source ./setup_lxplus.sh

arguments=(
  --job-id "${job_id}"
  --events "${events}"
  --nch-min "${nch_min}"
  --campaign-seed "${campaign_seed}"
  --output "${output_file}"
)
if [[ "${save_constituents}" == "1" ]]; then
  arguments+=(--save-constituents)
fi

./pp_run3_cp5_13p6 "${arguments[@]}"
