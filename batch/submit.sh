#!/usr/bin/env bash
# Submit pp_run3_cp5 condor jobs without editing the .sub file.
#
# Usage:
#   ./batch/submit.sh [options]
#
# Options:
#   --n-jobs          N      Number of jobs          (default: 10)
#   --events          N      Events per job           (default: 10000)
#   --nch-min         N      Min charged constituents (default: 60)
#   --campaign-seed   N      Campaign seed            (default: 246800099)
#   --job-flavour     STR    HTCondor job flavour     (default: workday)
#   -h, --help
#
# Examples:
#   # Test run (defaults):
#   ./batch/submit.sh
#
#   # Full high-Nch production:
#   ./batch/submit.sh --n-jobs 500 --events 50000 --campaign-seed 246800001
#
#   # Inclusive run:
#   ./batch/submit.sh --nch-min 0 --campaign-seed 136000001 --events 20000

set -euo pipefail

# --- Defaults ---
N_JOBS=10
EVENTS_PER_JOB=10000
NCH_MIN=60
CAMPAIGN_SEED=246800099
JOB_FLAVOUR=workday

# --- Parse arguments ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --n-jobs)          N_JOBS="$2";          shift 2 ;;
    --events)          EVENTS_PER_JOB="$2";  shift 2 ;;
    --nch-min)         NCH_MIN="$2";         shift 2 ;;
    --campaign-seed)   CAMPAIGN_SEED="$2";   shift 2 ;;
    --job-flavour)     JOB_FLAVOUR="$2";     shift 2 ;;
    -h|--help)
      sed -n '/^# Usage:/,/^[^#]/p' "$0" | grep '^#' | sed 's/^# \?//'
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# --- Ensure logs directory exists ---
mkdir -p logs

# --- Pre-create EOS output directory from login node (reliable EOS access) ---
EOS_DIR="/eos/cms/store/group/phys_heavyions/xiaoyul/13p6_pythia8/nchGE${NCH_MIN}_cs${CAMPAIGN_SEED}"
eos mkdir -p "${EOS_DIR}"
echo "  EOS dir ready: ${EOS_DIR}"

# --- Summary ---
echo "Submitting pp_run3_cp5 jobs:"
echo "  n_jobs          = ${N_JOBS}"
echo "  events_per_job  = ${EVENTS_PER_JOB}"
echo "  nch_min         = ${NCH_MIN}"
echo "  campaign_seed   = ${CAMPAIGN_SEED}"
echo "  job_flavour     = ${JOB_FLAVOUR}"
echo "  eos_dir         = /eos/cms/store/group/phys_heavyions/xiaoyul/13p6_pythia8/nchGE${NCH_MIN}_cs${CAMPAIGN_SEED}"
echo ""

condor_submit batch/submit_pp_run3_cp5.sub \
  n_jobs="${N_JOBS}"                        \
  events_per_job="${EVENTS_PER_JOB}"        \
  nch_min="${NCH_MIN}"                      \
  campaign_seed="${CAMPAIGN_SEED}"          \
  job_flavour="${JOB_FLAVOUR}"
