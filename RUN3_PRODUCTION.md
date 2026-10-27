# Standalone CMS Run-3 CP5 production at 13.6 TeV

`pp_run3_cp5_13p6.cc` generates standalone generator-level hard-QCD events
with the official CMS Run-3 CP5 settings for 13.6 TeV, FastJet anti-kT R=0.8
jets, and ROOT output. It does not require CMSSW.

## Physics definition

- pp collisions at sqrt(s) = 13.6 TeV.
- `HardQCD:all = on`, with default `pTHatMin = 500 GeV`.
- Jets are clustered from all stable visible final-state particles.
- Selected jets have `pT >= 550 GeV` and `|eta| <= 1.6`.
- Jet Nch counts charged constituents with `pT >= 0.3 GeV` and
  `|eta| <= 2.4`.
- `--nch-min 0` is inclusive in jet Nch.
- `--nch-min 60` retains events containing at least one selected jet with
  `Nch >= 60`. All selected jets in a retained event are written, preserving
  the original pT ordering. Use `genJetPassNchCut` to identify the jets that
  passed the skim.

The lightweight `trackTree` contains event metadata and per-jet vectors. Add
`--save-constituents` only when constituent-level output is required; it creates
a flat `constituentTree` keyed by `jobId`, `eventIndex`, and `jetRank`.

The one-entry `runInfo` tree stores the complete cuts, random seeds, generated
and saved counters, weight sums, and PYTHIA cross section. Those counters must
be retained when normalizing a filtered high-Nch sample.

## Random seeds

The PYTHIA seed is generated deterministically as

```
1 + ((campaignSeed + 48271 * jobId) mod 900000000)
```

The multiplier is coprime to PYTHIA's allowed seed range, so distinct Condor
job IDs cannot collide before the full range wraps. A retry of a job reproduces
the same events. This is preferable to timestamps, which can collide when many
jobs start together and make a production impossible to reproduce.

Keep `campaignSeed` fixed within one production. Choose and record a different
campaign seed for a statistically independent production. In particular, use a
different campaign seed if the inclusive and high-Nch productions should not
contain overlapping event streams.

## Build on lxplus

```bash
cd MC_Fragments
source batch/setup_lxplus.sh
make -f Makefile.run3
```

The setup defaults to an x86_64 EL9 LCG view. If needed, select another view:

```bash
export LCG_VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_XXX/x86_64-el9-gccXX-opt
source batch/setup_lxplus.sh
```

The environment must provide ROOT, FastJet, PYTHIA 8 built with LHAPDF6, and
the `NNPDF31_nnlo_as_0118` data set.

## Local validation commands

Inclusive test:

```bash
./pp_run3_cp5_13p6 \
  --job-id 0 \
  --events 1000 \
  --nch-min 0 \
  --campaign-seed 136000001 \
  --output test_inclusive.root
```

High-Nch test:

```bash
./pp_run3_cp5_13p6 \
  --job-id 0 \
  --events 10000 \
  --nch-min 60 \
  --campaign-seed 246800001 \
  --output test_nchGE60.root
```

## HTCondor on lxplus

Build the executable first, then create the return directories and submit:

```bash
mkdir -p logs outputs
condor_submit batch/submit_pp_run3_cp5.sub
```

The submit file defaults to ten jobs of ten thousand events with `nch_min = 0`.
Inspect the ROOT output and logs before increasing `n_jobs` or
`events_per_job`. For the later enriched production, set `nch_min = 60`, pick a
new recorded `campaign_seed`, and decide whether constituent output is needed.

Do not delete the `runInfo` trees when merging. Their generated-event counts,
weight sums, and cross sections are needed for an absolute comparison of the
13.0 and 13.6 TeV jet-Nch spectra.
