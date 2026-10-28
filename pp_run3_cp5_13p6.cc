// Standalone PYTHIA 8 + FastJet generator for CMS Run-3 CP5 pp events at
// sqrt(s) = 13.6 TeV.
//
// Tree format matches the Run-2 pp_highMultGen.cc output exactly:
//   - Event-level: all final-state particles (px, py, pz, m, pid, chg)
//   - Jet-level:   genJetPt/Eta/Phi/ChargedMultiplicity
//   - Constituent: genDau_pt/eta/phi/chg/pid (vector<vector<T>>)
//   - runInfo:     normalization and configuration metadata

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "TInterpreter.h"
#include "TNamed.h"
#include "TTree.h"
#include "TMath.h"

#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double     kEcmGeV         = 13600.0;
constexpr std::uint64_t kPythiaSeedRange = 900000000ULL;
constexpr std::uint64_t kSeedStride      =     48271ULL;

struct Options {
  std::uint64_t jobId        = 0;
  std::uint64_t campaignSeed = 136000001ULL;
  long long     events       = 10000;
  int           nchMin       = 0;
  double        pTHatMin     = 500.0;
  double        jetPtMin     = 550.0;
  double        jetAbsEtaMax = 1.6;
  double        constituentPtMin     = 0.3;
  double        constituentAbsEtaMax = 2.4;
  double        jetRadius    = 0.8;
  std::string   output;
};

void printUsage(const char* prog) {
  std::cout
    << "Usage: " << prog << " [options]\n\n"
    << "  --job-id N           Batch job identifier (default: 0)\n"
    << "  --campaign-seed N    Fixed seed identifying this production (default: 136000001)\n"
    << "  --events N           PYTHIA event attempts (default: 10000)\n"
    << "  --nch-min N          Min charged constituents in at least one jet; 0 = inclusive\n"
    << "  --pthat-min GEV      Hard-QCD pTHat minimum (default: 500)\n"
    << "  --jet-pt-min GEV     Jet pT minimum (default: 550)\n"
    << "  --jet-eta-max VALUE  Jet |eta| maximum (default: 1.6)\n"
    << "  --output FILE.root   Output file (default: auto-named)\n"
    << "  --help\n";
}

std::string requireValue(int argc, char* argv[], int& i) {
  if (i + 1 >= argc)
    throw std::runtime_error(std::string("Missing value after ") + argv[i]);
  return argv[++i];
}

Options parseOptions(int argc, char* argv[]) {
  Options o;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if      (a == "--help")          { printUsage(argv[0]); std::exit(0); }
    else if (a == "--job-id")        { o.jobId        = std::stoull(requireValue(argc,argv,i)); }
    else if (a == "--campaign-seed") { o.campaignSeed = std::stoull(requireValue(argc,argv,i)); }
    else if (a == "--events")        { o.events       = std::stoll (requireValue(argc,argv,i)); }
    else if (a == "--nch-min")       { o.nchMin       = std::stoi  (requireValue(argc,argv,i)); }
    else if (a == "--pthat-min")     { o.pTHatMin     = std::stod  (requireValue(argc,argv,i)); }
    else if (a == "--jet-pt-min")    { o.jetPtMin     = std::stod  (requireValue(argc,argv,i)); }
    else if (a == "--jet-eta-max")   { o.jetAbsEtaMax = std::stod  (requireValue(argc,argv,i)); }
    else if (a == "--output")        { o.output       = requireValue(argc,argv,i); }
    else throw std::runtime_error("Unknown argument: " + a);
  }
  if (o.output.empty()) {
    std::ostringstream s;
    s << "pp_run3_cp5_13p6_"
      << (o.nchMin == 0 ? "inclusive" : "nchGE" + std::to_string(o.nchMin))
      << "_job" << std::setw(6) << std::setfill('0') << o.jobId << ".root";
    o.output = s.str();
  }
  return o;
}

// Deterministic, collision-free seed across jobs:
//   seed = 1 + ((campaignSeed + 48271 * jobId) mod 900,000,000)
// 48271 is coprime to 900,000,000, so the mapping is bijective.
int makePythiaSeed(std::uint64_t campaignSeed, std::uint64_t jobId) {
  const std::uint64_t z =
    ((campaignSeed % kPythiaSeedRange) +
     (kSeedStride * (jobId % kPythiaSeedRange)) % kPythiaSeedRange)
    % kPythiaSeedRange;
  return static_cast<int>(z + 1ULL);
}

void configurePythia(Pythia8::Pythia& pythia, const Options& o, int seed) {
  const auto set = [&](const std::string& s) {
    if (!pythia.readString(s))
      throw std::runtime_error("PYTHIA rejected: " + s);
  };

  set("Beams:idA = 2212");
  set("Beams:idB = 2212");
  set("Beams:eCM = 13600.");

  // General CMS settings
  set("Check:epTolErr = 0.01");
  set("Beams:setProductionScalesFromLHEF = off");
  set("ParticleDecays:limitTau0 = on");
  set("ParticleDecays:tau0Max = 10");
  set("ParticleDecays:allowPhotonRadiation = on");

  // CMS CP5 tune — 13.6 TeV parameters.
  // PDF:pSet must come AFTER Tune:pp (tune selection can reset the PDF choice).
  set("Tune:pp = 14");
  set("Tune:ee = 7");
  set("MultipartonInteractions:ecmPow = 0.03344");
  set("MultipartonInteractions:bProfile = 2");
  set("MultipartonInteractions:pT0Ref = 1.41");
  set("MultipartonInteractions:coreRadius = 0.7634");
  set("MultipartonInteractions:coreFraction = 0.63");
  set("ColourReconnection:range = 5.176");
  set("SigmaTotal:zeroAXB = off");
  set("SpaceShower:alphaSorder = 2");
  set("SpaceShower:alphaSvalue = 0.118");
  set("SigmaProcess:alphaSvalue = 0.118");
  set("SigmaProcess:alphaSorder = 2");
  set("MultipartonInteractions:alphaSvalue = 0.118");
  set("MultipartonInteractions:alphaSorder = 2");
  set("TimeShower:alphaSorder = 2");
  set("TimeShower:alphaSvalue = 0.118");
  set("SigmaTotal:mode = 0");
  set("SigmaTotal:sigmaEl = 22.08");
  set("SigmaTotal:sigmaTot = 101.037");
  set("PDF:pSet = LHAPDF6:NNPDF31_nnlo_as_0118");

  // Hard process
  set("HardQCD:all = on");
  set("PhaseSpace:pTHatMin = " + std::to_string(o.pTHatMin));

  // Deterministic seed
  set("Random:setSeed = on");
  set("Random:seed = " + std::to_string(seed));
}

} // namespace

int main(int argc, char* argv[]) {

  // ROOT needs these dictionaries for vector<vector<T>> branches
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >",   "vector");

  Options options;
  try { options = parseOptions(argc, argv); }
  catch (const std::exception& e) {
    std::cerr << "Argument error: " << e.what() << "\n\n";
    printUsage(argv[0]); return 2;
  }

  const int pythiaSeed = makePythiaSeed(options.campaignSeed, options.jobId);
  std::cout << "Run-3 CP5 13.6 TeV production\n"
            << "  job ID:        " << options.jobId       << '\n'
            << "  campaign seed: " << options.campaignSeed << '\n'
            << "  PYTHIA seed:   " << pythiaSeed           << '\n'
            << "  events:        " << options.events       << '\n'
            << "  nch-min:       " << options.nchMin       << '\n'
            << "  output:        " << options.output       << std::endl;

  Pythia8::Pythia pythia;
  try { configurePythia(pythia, options, pythiaSeed); }
  catch (const std::exception& e) {
    std::cerr << "Config error: " << e.what() << std::endl; return 3;
  }
  if (!pythia.init()) {
    std::cerr << "PYTHIA init failed. Check LHAPDF6 and NNPDF31_nnlo_as_0118." << std::endl;
    return 3;
  }
  std::cout << "Effective PDF:pSet = " << pythia.settings.word("PDF:pSet")
            << "\nEffective MPI ecmPow = "
            << pythia.settings.parm("MultipartonInteractions:ecmPow") << std::endl;

  // ── Output file ──────────────────────────────────────────────────────────
  TFile* f = TFile::Open(options.output.c_str(), "RECREATE");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot create " << options.output << std::endl; return 4;
  }

  // ── Tree branches — identical to Run-2 pp_highMultGen.cc ─────────────────
  float NchHadrons = 0;
  float weight     = 0;
  int   processCode = 0;
  int   nMPI        = 0;
  // All final-state particles
  std::vector<float> px, py, pz, m;
  std::vector<int>   pid, chg;
  // DIS placeholders (kept for branch-name compatibility; always 0 in pp)
  float Q2 = 0, W2 = 0, x = 0, y = 0;
  // Jet kinematics
  std::vector<float> genJetPt, genJetEta, genJetPhi, genJetMass;
  std::vector<int>   genJetChargedMultiplicity;
  // Jet constituents (vector<vector<T>>)
  std::vector<std::vector<float>> genDau_pt, genDau_eta, genDau_phi;
  std::vector<std::vector<int>>   genDau_chg, genDau_pid;

  TTree* trackTree = new TTree("trackTree", "Run-3 CP5 generator-level events");
  trackTree->Branch("NchHadrons",              &NchHadrons);
  trackTree->Branch("weight",                  &weight);
  trackTree->Branch("processCode",             &processCode);
  trackTree->Branch("nMPI",                    &nMPI);
  trackTree->Branch("px",                      &px);
  trackTree->Branch("py",                      &py);
  trackTree->Branch("pz",                      &pz);
  trackTree->Branch("m",                       &m);
  trackTree->Branch("pid",                     &pid);
  trackTree->Branch("chg",                     &chg);
  trackTree->Branch("Q2",                      &Q2);
  trackTree->Branch("W2",                      &W2);
  trackTree->Branch("x",                       &x);
  trackTree->Branch("y",                       &y);
  trackTree->Branch("genJetPt",                &genJetPt);
  trackTree->Branch("genJetEta",               &genJetEta);
  trackTree->Branch("genJetPhi",               &genJetPhi);
  trackTree->Branch("genJetMass",               &genJetMass);
  trackTree->Branch("genJetChargedMultiplicity",&genJetChargedMultiplicity);
  trackTree->Branch("genDau_pt",               &genDau_pt);
  trackTree->Branch("genDau_eta",              &genDau_eta);
  trackTree->Branch("genDau_phi",              &genDau_phi);
  trackTree->Branch("genDau_chg",              &genDau_chg);
  trackTree->Branch("genDau_pid",              &genDau_pid);

  // ── FastJet ──────────────────────────────────────────────────────────────
  const fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, options.jetRadius);

  // ── Counters for runInfo ─────────────────────────────────────────────────
  long long generatedEvents        = 0;
  long long savedEvents            = 0;
  long long selectedJets           = 0;
  long long jetsPassingNch         = 0;
  double    sumGeneratedWeights    = 0.0;
  double    sumSavedWeights        = 0.0;
  double    sumSelectedJetWeights  = 0.0;
  double    sumPassingJetWeights   = 0.0;

  // ── Event loop ───────────────────────────────────────────────────────────
  for (long long attempt = 0; attempt < options.events; ++attempt) {
    if (attempt % 10000 == 0)
      std::cout << "Event attempt " << attempt << " / " << options.events << std::endl;
    if (!pythia.next()) continue;

    ++generatedEvents;
    const double evtWeight = pythia.info.weight();
    sumGeneratedWeights += evtWeight;

    // ── Collect all final-state particles ─────────────────────────────────
    std::vector<fastjet::PseudoJet> particles;
    float nchHadronsEv = 0;
    for (int i = 0; i < pythia.event.size(); ++i) {
      const Pythia8::Particle& p = pythia.event[i];
      if (!p.isFinal()) continue;
      if (p.isCharged() && p.isHadron()) ++nchHadronsEv;
      // All final-state for event branches
      px.push_back(p.px()); py.push_back(p.py());
      pz.push_back(p.pz()); m.push_back(p.m());
      pid.push_back(p.id()); chg.push_back(static_cast<int>(p.charge()));
      // For jet clustering: visible particles with basic kinematic cuts
      if (!p.isVisible()) continue;
      if (TMath::Abs(p.eta()) > options.constituentAbsEtaMax || p.pT() < options.constituentPtMin) continue;
      fastjet::PseudoJet pj(p.px(), p.py(), p.pz(), p.e());
      pj.set_user_index(i);
      particles.push_back(pj);
    }

    // ── Cluster jets ──────────────────────────────────────────────────────
    if (particles.empty()) {
      px.clear(); py.clear(); pz.clear(); m.clear(); pid.clear(); chg.clear();
      continue;
    }
    const fastjet::ClusterSequence cs(particles, jetDef);
    const std::vector<fastjet::PseudoJet> allJets =
      fastjet::sorted_by_pt(cs.inclusive_jets(options.jetPtMin));

    // Build analysis jets with Nch
    std::vector<fastjet::PseudoJet> analysisJets;
    std::vector<int>                analysisNch;
    std::vector<std::vector<float>> dauPt, dauEta, dauPhi;
    std::vector<std::vector<int>>   dauChg, dauPid;

    for (const auto& jet : allJets) {
      if (TMath::Abs(jet.eta()) > options.jetAbsEtaMax) continue;
      int chMult = 0;
      std::vector<float> tPt, tEta, tPhi;
      std::vector<int>   tChg, tPid;
      for (const auto& c : jet.constituents()) {
        const Pythia8::Particle& cp = pythia.event[c.user_index()];
        if (cp.isCharged()) ++chMult;
        tPt.push_back(c.pt());
        tEta.push_back(c.eta());
        tPhi.push_back(c.phi());
        tChg.push_back(static_cast<int>(cp.charge()));
        tPid.push_back(cp.id());
      }
      analysisJets.push_back(jet);
      analysisNch.push_back(chMult);
      dauPt.push_back(tPt); dauEta.push_back(tEta); dauPhi.push_back(tPhi);
      dauChg.push_back(tChg); dauPid.push_back(tPid);
    }

    if (analysisJets.empty()) {
      px.clear(); py.clear(); pz.clear(); m.clear(); pid.clear(); chg.clear();
      continue;
    }

    // Apply Nch filter: require at least one jet with Nch >= nchMin
    int nPassingJets = 0;
    for (int nch : analysisNch) if (nch >= options.nchMin) ++nPassingJets;
    if (options.nchMin > 0 && nPassingJets == 0) {
      px.clear(); py.clear(); pz.clear(); m.clear(); pid.clear(); chg.clear();
      continue;
    }

    // ── Fill tree ─────────────────────────────────────────────────────────
    NchHadrons  = nchHadronsEv;
    weight      = static_cast<float>(evtWeight);
    processCode = pythia.info.code();
    nMPI        = pythia.info.nMPI();
    Q2 = 0; W2 = 0; x = 0; y = 0; // not applicable in pp

    for (std::size_t j = 0; j < analysisJets.size(); ++j) {
      genJetPt.push_back(static_cast<float>(analysisJets[j].pt()));
      genJetEta.push_back(static_cast<float>(analysisJets[j].eta()));
      genJetPhi.push_back(static_cast<float>(analysisJets[j].phi()));
      genJetMass.push_back(static_cast<float>(analysisJets[j].m()));
      genJetChargedMultiplicity.push_back(analysisNch[j]);
      genDau_pt.push_back(dauPt[j]);
      genDau_eta.push_back(dauEta[j]);
      genDau_phi.push_back(dauPhi[j]);
      genDau_chg.push_back(dauChg[j]);
      genDau_pid.push_back(dauPid[j]);
    }

    trackTree->Fill();

    ++savedEvents;
    selectedJets        += static_cast<long long>(analysisJets.size());
    jetsPassingNch      += nPassingJets;
    sumSavedWeights     += evtWeight;
    sumSelectedJetWeights += evtWeight * analysisJets.size();
    sumPassingJetWeights  += evtWeight * nPassingJets;

    // Clear per-event vectors
    NchHadrons = 0; weight = 0; processCode = 0; nMPI = 0;
    px.clear(); py.clear(); pz.clear(); m.clear(); pid.clear(); chg.clear();
    genJetPt.clear(); genJetEta.clear(); genJetPhi.clear(); genJetMass.clear();
    genJetChargedMultiplicity.clear();
    genDau_pt.clear(); genDau_eta.clear(); genDau_phi.clear();
    genDau_chg.clear(); genDau_pid.clear();
  }

  pythia.stat();

  // ── runInfo tree ─────────────────────────────────────────────────────────
  TTree runInfo("runInfo", "Production configuration and normalization");
  double    ecmGeV_val               = kEcmGeV;
  auto      jobId_val                = options.jobId;
  auto      campaignSeed_val         = options.campaignSeed;
  int       pythiaSeed_val           = pythiaSeed;
  long long attemptedEvents_val      = options.events;
  int       nchMin_val               = options.nchMin;
  double    pTHatMin_val             = options.pTHatMin;
  double    jetPtMin_val             = options.jetPtMin;
  double    jetAbsEtaMax_val         = options.jetAbsEtaMax;
  double    constituentPtMin_val     = options.constituentPtMin;
  double    constituentAbsEtaMax_val = options.constituentAbsEtaMax;
  double    jetRadius_val            = options.jetRadius;
  double    sigmaGenMb               = pythia.info.sigmaGen();
  double    sigmaErrMb               = pythia.info.sigmaErr();

  runInfo.Branch("ecmGeV",                &ecmGeV_val);
  runInfo.Branch("jobId",                 &jobId_val);
  runInfo.Branch("campaignSeed",          &campaignSeed_val);
  runInfo.Branch("pythiaSeed",            &pythiaSeed_val);
  runInfo.Branch("attemptedEvents",       &attemptedEvents_val);
  runInfo.Branch("generatedEvents",       &generatedEvents);
  runInfo.Branch("savedEvents",           &savedEvents);
  runInfo.Branch("selectedJets",          &selectedJets);
  runInfo.Branch("jetsPassingNch",        &jetsPassingNch);
  runInfo.Branch("sumGeneratedWeights",   &sumGeneratedWeights);
  runInfo.Branch("sumSavedWeights",       &sumSavedWeights);
  runInfo.Branch("sumSelectedJetWeights", &sumSelectedJetWeights);
  runInfo.Branch("sumPassingJetWeights",  &sumPassingJetWeights);
  runInfo.Branch("sigmaGenMb",            &sigmaGenMb);
  runInfo.Branch("sigmaErrMb",            &sigmaErrMb);
  runInfo.Branch("nchMin",               &nchMin_val);
  runInfo.Branch("pTHatMin",             &pTHatMin_val);
  runInfo.Branch("jetPtMin",             &jetPtMin_val);
  runInfo.Branch("jetAbsEtaMax",         &jetAbsEtaMax_val);
  runInfo.Branch("constituentPtMin",     &constituentPtMin_val);
  runInfo.Branch("constituentAbsEtaMax", &constituentAbsEtaMax_val);
  runInfo.Branch("jetRadius",            &jetRadius_val);
  runInfo.Fill();

  // ── Write and close ───────────────────────────────────────────────────────
  TNamed tuneMeta("tune",   "CMS CP5 Run3 ECM 13.6 TeV");
  TNamed pdfMeta ("pdfSet", pythia.settings.word("PDF:pSet").c_str());
  f->cd();
  tuneMeta.Write(); pdfMeta.Write();
  trackTree->Write();
  runInfo.Write();
  f->Close();

  std::cout << "Generated " << generatedEvents << " events; saved " << savedEvents
            << " events with " << selectedJets << " jets ("
            << jetsPassingNch << " pass Nch >= " << options.nchMin << ")." << std::endl;
  return 0;
}
