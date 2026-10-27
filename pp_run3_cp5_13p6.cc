// Standalone PYTHIA 8 + FastJet generator for CMS Run-3 CP5 pp events at
// sqrt(s) = 13.6 TeV.  The output is generator-level ROOT only.

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "TNamed.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kEcmGeV = 13600.0;
constexpr std::uint64_t kPythiaSeedRange = 900000000ULL;
constexpr std::uint64_t kSeedStride = 48271ULL;

struct Options {
  std::uint64_t jobId = 0;
  std::uint64_t campaignSeed = 136000001ULL;
  long long events = 10000;
  int nchMin = 0;
  double pTHatMin = 500.0;
  double jetPtMin = 550.0;
  double jetAbsEtaMax = 1.6;
  double constituentPtMin = 0.3;
  double constituentAbsEtaMax = 2.4;
  double jetRadius = 0.8;
  bool saveConstituents = false;
  std::string output;
};

void printUsage(const char* program) {
  std::cout
      << "Usage: " << program << " [options]\n\n"
      << "  --job-id N                 Batch job identifier (default: 0)\n"
      << "  --campaign-seed N          Fixed seed identifying this production\n"
      << "                             (default: 136000001)\n"
      << "  --events N                 PYTHIA event attempts (default: 10000)\n"
      << "  --nch-min N                Require at least one selected jet with\n"
      << "                             Nch >= N; 0 is inclusive (default: 0)\n"
      << "  --pthat-min GEV            Hard-QCD pTHat minimum (default: 500)\n"
      << "  --jet-pt-min GEV           Analysis-jet pT minimum (default: 550)\n"
      << "  --jet-eta-max VALUE        Analysis-jet |eta| maximum (default: 1.6)\n"
      << "  --save-constituents        Write the optional constituentTree\n"
      << "  --output FILE.root         Output path (otherwise generated)\n"
      << "  --help                     Show this message\n";
}

std::string requireValue(int argc, char* argv[], int& index) {
  if (index + 1 >= argc) {
    throw std::runtime_error(std::string("Missing value after ") + argv[index]);
  }
  return argv[++index];
}

Options parseOptions(int argc, char* argv[]) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--help") {
      printUsage(argv[0]);
      std::exit(0);
    } else if (argument == "--job-id") {
      options.jobId = std::stoull(requireValue(argc, argv, i));
    } else if (argument == "--campaign-seed") {
      options.campaignSeed = std::stoull(requireValue(argc, argv, i));
    } else if (argument == "--events") {
      options.events = std::stoll(requireValue(argc, argv, i));
    } else if (argument == "--nch-min") {
      options.nchMin = std::stoi(requireValue(argc, argv, i));
    } else if (argument == "--pthat-min") {
      options.pTHatMin = std::stod(requireValue(argc, argv, i));
    } else if (argument == "--jet-pt-min") {
      options.jetPtMin = std::stod(requireValue(argc, argv, i));
    } else if (argument == "--jet-eta-max") {
      options.jetAbsEtaMax = std::stod(requireValue(argc, argv, i));
    } else if (argument == "--save-constituents") {
      options.saveConstituents = true;
    } else if (argument == "--output") {
      options.output = requireValue(argc, argv, i);
    } else {
      throw std::runtime_error("Unknown argument: " + argument);
    }
  }

  if (options.events <= 0) throw std::runtime_error("--events must be positive");
  if (options.nchMin < 0) throw std::runtime_error("--nch-min cannot be negative");
  if (options.pTHatMin <= 0.0) throw std::runtime_error("--pthat-min must be positive");
  if (options.jetPtMin <= 0.0) throw std::runtime_error("--jet-pt-min must be positive");
  if (options.jetAbsEtaMax <= 0.0) throw std::runtime_error("--jet-eta-max must be positive");

  if (options.output.empty()) {
    std::ostringstream name;
    name << "pp_run3_cp5_13p6_";
    if (options.nchMin == 0) {
      name << "inclusive";
    } else {
      name << "nchGE" << options.nchMin;
    }
    name << "_job" << std::setw(6) << std::setfill('0') << options.jobId << ".root";
    options.output = name.str();
  }
  return options;
}

// 48271 is coprime to 900,000,000, so this affine mapping is one-to-one
// over a full PYTHIA seed period.  Different job IDs therefore cannot collide
// until the job ID itself wraps modulo the complete allowed seed range.
int makePythiaSeed(std::uint64_t campaignSeed, std::uint64_t jobId) {
  const std::uint64_t zeroBased =
      ((campaignSeed % kPythiaSeedRange) +
       (kSeedStride * (jobId % kPythiaSeedRange)) % kPythiaSeedRange) %
      kPythiaSeedRange;
  return static_cast<int>(zeroBased + 1ULL);
}

void configurePythia(Pythia8::Pythia& pythia, const Options& options,
                     int pythiaSeed) {
  const auto set = [&pythia](const std::string& setting) {
    if (!pythia.readString(setting)) {
      throw std::runtime_error("PYTHIA rejected setting: " + setting);
    }
  };

  set("Beams:idA = 2212");
  set("Beams:idB = 2212");
  set("Beams:eCM = 13600.");

  // General CMS generator settings.
  set("Check:epTolErr = 0.01");
  set("Beams:setProductionScalesFromLHEF = off");
  set("ParticleDecays:limitTau0 = on");
  set("ParticleDecays:tau0Max = 10");
  set("ParticleDecays:allowPhotonRadiation = on");

  // CMS CP5, using the Run-3 13.6 TeV total-cross-section parameters.
  // Keep PDF:pSet after Tune:pp because selecting a pp tune can reset the PDF.
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

  // Process settings come after the tune, matching the ordering used by CMS.
  set("HardQCD:all = on");
  set("PhaseSpace:pTHatMin = " + std::to_string(options.pTHatMin));

  set("Random:setSeed = on");
  set("Random:seed = " + std::to_string(pythiaSeed));
}

}  // namespace

int main(int argc, char* argv[]) {
  Options options;
  try {
    options = parseOptions(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "Argument error: " << error.what() << "\n\n";
    printUsage(argv[0]);
    return 2;
  }

  const int pythiaSeed = makePythiaSeed(options.campaignSeed, options.jobId);
  std::cout << "Run-3 CP5 13.6 TeV production configuration\n"
            << "  job ID:          " << options.jobId << '\n'
            << "  campaign seed:   " << options.campaignSeed << '\n'
            << "  PYTHIA seed:     " << pythiaSeed << '\n'
            << "  attempted events:" << options.events << '\n'
            << "  Nch minimum:     " << options.nchMin << '\n'
            << "  output:          " << options.output << std::endl;

  Pythia8::Pythia pythia;
  try {
    configurePythia(pythia, options, pythiaSeed);
  } catch (const std::exception& error) {
    std::cerr << "Configuration error: " << error.what() << std::endl;
    return 3;
  }
  if (!pythia.init()) {
    std::cerr << "PYTHIA initialization failed. Check LHAPDF6 and the "
                 "NNPDF31_nnlo_as_0118 data set."
              << std::endl;
    return 3;
  }
  std::cout << "Effective PDF:pSet = " << pythia.settings.word("PDF:pSet")
            << "\nEffective MPI ecmPow = "
            << pythia.settings.parm("MultipartonInteractions:ecmPow")
            << std::endl;

  TFile outputFile(options.output.c_str(), "RECREATE");
  if (outputFile.IsZombie()) {
    std::cerr << "Cannot create output file " << options.output << std::endl;
    return 4;
  }

  // One entry per saved event. Jets remain sorted by pT, and all selected jets
  // are retained even for a high-Nch skim.
  TTree trackTree("trackTree", "Run-3 CP5 generator-level jet events");
  unsigned long long treeJobId = options.jobId;
  long long eventIndex = 0;
  double weight = 0.0;
  double pTHat = 0.0;
  int processCode = 0;
  int nMPI = 0;
  int eventNch = 0;
  int nJetsPassingNch = 0;
  std::vector<float> genJetPt;
  std::vector<float> genJetEta;
  std::vector<float> genJetPhi;
  std::vector<float> genJetMass;
  std::vector<int> genJetChargedMultiplicity;
  std::vector<int> genJetPassNchCut;

  trackTree.Branch("jobId", &treeJobId);
  trackTree.Branch("eventIndex", &eventIndex);
  trackTree.Branch("weight", &weight);
  trackTree.Branch("pTHat", &pTHat);
  trackTree.Branch("processCode", &processCode);
  trackTree.Branch("nMPI", &nMPI);
  trackTree.Branch("eventNch", &eventNch);
  trackTree.Branch("nJetsPassingNch", &nJetsPassingNch);
  trackTree.Branch("genJetPt", &genJetPt);
  trackTree.Branch("genJetEta", &genJetEta);
  trackTree.Branch("genJetPhi", &genJetPhi);
  trackTree.Branch("genJetMass", &genJetMass);
  trackTree.Branch("genJetChargedMultiplicity", &genJetChargedMultiplicity);
  trackTree.Branch("genJetPassNchCut", &genJetPassNchCut);

  // Optional flat constituent tree. Flat rows avoid nested-vector dictionary
  // issues and can be joined to trackTree with jobId, eventIndex, and jetRank.
  TTree constituentTree("constituentTree", "Selected generator-jet constituents");
  unsigned long long constituentJobId = options.jobId;
  long long constituentEventIndex = 0;
  int constituentJetRank = 0;
  int constituentPid = 0;
  float constituentCharge = 0.0F;
  float constituentPt = 0.0F;
  float constituentEta = 0.0F;
  float constituentPhi = 0.0F;
  float constituentMass = 0.0F;
  bool constituentCountsTowardNch = false;
  if (options.saveConstituents) {
    constituentTree.Branch("jobId", &constituentJobId);
    constituentTree.Branch("eventIndex", &constituentEventIndex);
    constituentTree.Branch("jetRank", &constituentJetRank);
    constituentTree.Branch("pid", &constituentPid);
    constituentTree.Branch("charge", &constituentCharge);
    constituentTree.Branch("pt", &constituentPt);
    constituentTree.Branch("eta", &constituentEta);
    constituentTree.Branch("phi", &constituentPhi);
    constituentTree.Branch("mass", &constituentMass);
    constituentTree.Branch("countsTowardNch", &constituentCountsTowardNch);
  }

  const fastjet::JetDefinition jetDefinition(fastjet::antikt_algorithm,
                                              options.jetRadius);
  long long generatedEvents = 0;
  long long savedEvents = 0;
  long long selectedJets = 0;
  long long jetsPassingNch = 0;
  double sumGeneratedWeights = 0.0;
  double sumSavedWeights = 0.0;
  double sumSelectedJetWeights = 0.0;
  double sumPassingJetWeights = 0.0;

  for (long long attempt = 0; attempt < options.events; ++attempt) {
    if (attempt % 10000 == 0) {
      std::cout << "Event attempt " << attempt << " / " << options.events
                << std::endl;
    }
    if (!pythia.next()) continue;

    ++generatedEvents;
    const double eventWeight = pythia.info.weight();
    sumGeneratedWeights += eventWeight;

    std::vector<fastjet::PseudoJet> particles;
    particles.reserve(pythia.event.size());
    int currentEventNch = 0;
    for (int index = 1; index < pythia.event.size(); ++index) {
      const Pythia8::Particle& particle = pythia.event[index];
      if (!particle.isFinal() || !particle.isVisible()) continue;

      fastjet::PseudoJet pseudoJet(particle.px(), particle.py(), particle.pz(),
                                   particle.e());
      pseudoJet.set_user_index(index);
      particles.push_back(pseudoJet);

      if (particle.isCharged() && particle.pT() >= options.constituentPtMin &&
          std::abs(particle.eta()) <= options.constituentAbsEtaMax) {
        ++currentEventNch;
      }
    }
    if (particles.empty()) continue;

    const fastjet::ClusterSequence clusterSequence(particles, jetDefinition);
    const std::vector<fastjet::PseudoJet> allJets =
        fastjet::sorted_by_pt(clusterSequence.inclusive_jets(options.jetPtMin));

    std::vector<fastjet::PseudoJet> analysisJets;
    std::vector<int> analysisJetNch;
    for (const fastjet::PseudoJet& jet : allJets) {
      if (std::abs(jet.eta()) > options.jetAbsEtaMax) continue;
      int jetNch = 0;
      for (const fastjet::PseudoJet& constituent : jet.constituents()) {
        const Pythia8::Particle& particle =
            pythia.event[constituent.user_index()];
        if (particle.isCharged() && constituent.pt() >= options.constituentPtMin &&
            std::abs(constituent.eta()) <= options.constituentAbsEtaMax) {
          ++jetNch;
        }
      }
      analysisJets.push_back(jet);
      analysisJetNch.push_back(jetNch);
    }
    if (analysisJets.empty()) continue;

    const int currentPassingJets = static_cast<int>(std::count_if(
        analysisJetNch.begin(), analysisJetNch.end(),
        [&options](int nch) { return nch >= options.nchMin; }));
    if (options.nchMin > 0 && currentPassingJets == 0) continue;

    eventIndex = attempt;
    weight = eventWeight;
    pTHat = pythia.info.pTHat();
    processCode = pythia.info.code();
    nMPI = pythia.info.nMPI();
    eventNch = currentEventNch;
    nJetsPassingNch = currentPassingJets;
    genJetPt.clear();
    genJetEta.clear();
    genJetPhi.clear();
    genJetMass.clear();
    genJetChargedMultiplicity.clear();
    genJetPassNchCut.clear();

    for (std::size_t jetIndex = 0; jetIndex < analysisJets.size(); ++jetIndex) {
      const fastjet::PseudoJet& jet = analysisJets[jetIndex];
      const int jetNch = analysisJetNch[jetIndex];
      genJetPt.push_back(static_cast<float>(jet.pt()));
      genJetEta.push_back(static_cast<float>(jet.eta()));
      genJetPhi.push_back(static_cast<float>(jet.phi_std()));
      genJetMass.push_back(static_cast<float>(jet.m()));
      genJetChargedMultiplicity.push_back(jetNch);
      genJetPassNchCut.push_back(jetNch >= options.nchMin ? 1 : 0);

      if (!options.saveConstituents) continue;
      constituentEventIndex = attempt;
      constituentJetRank = static_cast<int>(jetIndex);
      for (const fastjet::PseudoJet& constituent : jet.constituents()) {
        const Pythia8::Particle& particle =
            pythia.event[constituent.user_index()];
        constituentPid = particle.id();
        constituentCharge = static_cast<float>(particle.charge());
        constituentPt = static_cast<float>(constituent.pt());
        constituentEta = static_cast<float>(constituent.eta());
        constituentPhi = static_cast<float>(constituent.phi_std());
        constituentMass = static_cast<float>(constituent.m());
        constituentCountsTowardNch =
            particle.isCharged() && constituent.pt() >= options.constituentPtMin &&
            std::abs(constituent.eta()) <= options.constituentAbsEtaMax;
        constituentTree.Fill();
      }
    }

    trackTree.Fill();
    ++savedEvents;
    selectedJets += static_cast<long long>(analysisJets.size());
    jetsPassingNch += currentPassingJets;
    sumSavedWeights += eventWeight;
    sumSelectedJetWeights += eventWeight * static_cast<double>(analysisJets.size());
    sumPassingJetWeights += eventWeight * static_cast<double>(currentPassingJets);
  }

  pythia.stat();

  TTree runInfo("runInfo", "Production configuration and normalization");
  double ecmGeV = kEcmGeV;
  unsigned long long jobId = options.jobId;
  unsigned long long campaignSeed = options.campaignSeed;
  int storedPythiaSeed = pythiaSeed;
  long long attemptedEvents = options.events;
  int nchMin = options.nchMin;
  double pTHatMin = options.pTHatMin;
  double jetPtMin = options.jetPtMin;
  double jetAbsEtaMax = options.jetAbsEtaMax;
  double constituentPtMin = options.constituentPtMin;
  double constituentAbsEtaMax = options.constituentAbsEtaMax;
  double jetRadius = options.jetRadius;
  bool savedConstituents = options.saveConstituents;
  double sigmaGenMb = pythia.info.sigmaGen();
  double sigmaErrMb = pythia.info.sigmaErr();

  runInfo.Branch("ecmGeV", &ecmGeV);
  runInfo.Branch("jobId", &jobId);
  runInfo.Branch("campaignSeed", &campaignSeed);
  runInfo.Branch("pythiaSeed", &storedPythiaSeed);
  runInfo.Branch("attemptedEvents", &attemptedEvents);
  runInfo.Branch("generatedEvents", &generatedEvents);
  runInfo.Branch("savedEvents", &savedEvents);
  runInfo.Branch("selectedJets", &selectedJets);
  runInfo.Branch("jetsPassingNch", &jetsPassingNch);
  runInfo.Branch("sumGeneratedWeights", &sumGeneratedWeights);
  runInfo.Branch("sumSavedWeights", &sumSavedWeights);
  runInfo.Branch("sumSelectedJetWeights", &sumSelectedJetWeights);
  runInfo.Branch("sumPassingJetWeights", &sumPassingJetWeights);
  runInfo.Branch("sigmaGenMb", &sigmaGenMb);
  runInfo.Branch("sigmaErrMb", &sigmaErrMb);
  runInfo.Branch("nchMin", &nchMin);
  runInfo.Branch("pTHatMin", &pTHatMin);
  runInfo.Branch("jetPtMin", &jetPtMin);
  runInfo.Branch("jetAbsEtaMax", &jetAbsEtaMax);
  runInfo.Branch("constituentPtMin", &constituentPtMin);
  runInfo.Branch("constituentAbsEtaMax", &constituentAbsEtaMax);
  runInfo.Branch("jetRadius", &jetRadius);
  runInfo.Branch("savedConstituents", &savedConstituents);
  runInfo.Fill();

  outputFile.cd();
  TNamed tuneMetadata("tune", "CMS CP5 Run3 ECM 13.6 TeV");
  TNamed pdfMetadata("pdfSet", pythia.settings.word("PDF:pSet").c_str());
  tuneMetadata.Write();
  pdfMetadata.Write();
  trackTree.Write();
  runInfo.Write();
  if (options.saveConstituents) constituentTree.Write();
  outputFile.Close();

  std::cout << "Generated " << generatedEvents << " events; saved " << savedEvents
            << " events containing " << selectedJets << " selected jets ("
            << jetsPassingNch << " pass Nch >= " << options.nchMin << ")."
            << std::endl;
  return 0;
}
