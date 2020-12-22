#include "AliAnalysisTaskDHFeCorr.h"


#include <stdexcept>

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
//#include "DHFe/DHFeConfig.h"
//#include "DHFe/DHFeNonHFe.h"
//#include "DHFe/DHFeOrigin.h"
#include "TChain.h"
#include "TFile.h"
#include "TGrid.h"

/*
namespace mdl = dhfe::model;
namespace qa = dhfe::qa;
namespace nhfe = dhfe::non_hfe;
namespace org = dhfe::origin;
namespace yaml = dhfe::yaml;
namespace cfg = dhfe::config;
namespace sel = dhfe::selection;


dhfe::config::DHFeTaskConfig::DHFeTaskConfig(
    const PWG::Tools::AliYAMLConfiguration &config) {
  const std::string top_level = "task";

  config.GetProperty(std::vector<std::string>({top_level, "mc_mode"}), fIsMC,
                     true);

  config.GetProperty(std::vector<std::string>({top_level, "qa_histograms"}),
                     fSaveHistograms, true);

  config.GetProperty(std::vector<std::string>({top_level, "only_efficiency"}),
                     fIsEffMode, true);

  config.GetProperty(std::vector<std::string>({top_level, "process_electron"}),
                     fProcessElectron, true);

  config.GetProperty(std::vector<std::string>({top_level, "process_dmeson"}),
                     fProcessDMeson, true);

  config.GetProperty(std::vector<std::string>({top_level, "save_event"}),
                     fSaveEvent, true);
}
*/

ClassImp(AliAnalysisTaskDHFeCorr);

AliAnalysisTaskDHFeCorr::AliAnalysisTaskDHFeCorr() : AliAnalysisTaskSE() {;}

AliAnalysisTaskDHFeCorr::~AliAnalysisTaskDHFeCorr() {;} 


AliAnalysisTaskDHFeCorr::AliAnalysisTaskDHFeCorr(const char *name,
                                                 unsigned int trigger,
                                                 const char *config,
                                                 const char *default_config)
    : AliAnalysisTaskSE(name) , 
    fDefaultConfig{std::string(default_config)},
    fUserConfig{std::string(config)} 
    {

  DefineIO();

  SelectCollisionCandidates(trigger);
  //ConfigureFromYaml();
  ConnectToAnalysisManager();
  CreateAndConnectIO(name);
}


void AliAnalysisTaskDHFeCorr::ConnectToAnalysisManager() {
  auto *manager = AliAnalysisManager::GetAnalysisManager();

  if (!manager)
    throw std::runtime_error("Not connected to the AnalysisManager.");

  if (!manager->GetInputEventHandler())
    throw std::runtime_error("Not connect to the InputEventHandler.");

  manager->AddTask(this);
}

void AliAnalysisTaskDHFeCorr::DefineIO() {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
  DefineOutput(7, TTree::Class());
  DefineOutput(8, TTree::Class());
}


void AliAnalysisTaskDHFeCorr::CreateAndConnectIO(const std::string &name) {
  auto *manager = AliAnalysisManager::GetAnalysisManager();

  if (!manager)
    throw std::runtime_error("Not connected to the AnalysisManager.");

  auto input = manager->GetCommonInputContainer();
  if (!input) throw std::runtime_error("No InputContainer.");

  this->ConnectInput(0, input);

  std::string file_output = AliAnalysisManager::GetCommonFileName();
  file_output += ":DHFeCorrelation_" + name;

  AliAnalysisDataContainer *electron_container = manager->CreateContainer(
      "ElectronTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
      file_output.c_str());
  this->ConnectOutput(1, electron_container);

  AliAnalysisDataContainer *dmeson_container = manager->CreateContainer(
      "DMesonTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
      file_output.c_str());
  this->ConnectOutput(2, dmeson_container);

  AliAnalysisDataContainer *event_container = manager->CreateContainer(
      "EventTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
      file_output.c_str());
  this->ConnectOutput(3, event_container);

  AliAnalysisDataContainer *event_qa_container = manager->CreateContainer(
      "EventsQA", TList::Class(), AliAnalysisManager::kOutputContainer,
      file_output.c_str());
  this->ConnectOutput(4, event_qa_container);

  AliAnalysisDataContainer *electron_qa_container = manager->CreateContainer(
      "ElectronQA", TList::Class(), AliAnalysisManager::kOutputContainer,
      file_output.c_str());
  this->ConnectOutput(5, electron_qa_container);

  AliAnalysisDataContainer *dmeson_qa_container = manager->CreateContainer(
      "DMesonQA", TList::Class(), AliAnalysisManager::kOutputContainer,
      file_output.c_str());
  this->ConnectOutput(6, dmeson_qa_container);

  AliAnalysisDataContainer *electron_mc_container = manager->CreateContainer(
      "ElectronTreeMC", TTree::Class(), AliAnalysisManager::kOutputContainer,
      file_output.c_str());
  this->ConnectOutput(7, electron_mc_container);

  AliAnalysisDataContainer *dmeson_mc_container = manager->CreateContainer(
      "DMesonTreeMC", TTree::Class(), AliAnalysisManager::kOutputContainer,
      file_output.c_str());
  this->ConnectOutput(8, dmeson_mc_container);
}

/*
void AliAnalysisTaskDHFeCorr::ConfigureFromYaml() {
  TGrid::Connect("alien://");

  fYAMLConfig.AddConfiguration(fDefaultConfig.c_str(), "default_configuration");
  fYAMLConfig.AddConfiguration(fUserConfig, std::string(fName));

  if (!fYAMLConfig.Initialize())
    throw std::runtime_error("Not possible to prepare the YAML configuration."); 

  fConfig = cfg::DHFeTaskConfig(fYAMLConfig);
  if (fConfig.ProcessElectron()) {
    fMainESelection = sel::ElectronSelection("main_electron", fYAMLConfig);
    std::cout << fMainESelection << std::endl;

    fPartnerESelection =
        sel::ElectronSelection("partner_electron", fYAMLConfig);
    std::cout << fPartnerESelection << std::endl;

    fMainEQAConfig = qa::ElectronQAConfig("main_electron", fYAMLConfig);
    fPartnerEQAConfig = qa::ElectronQAConfig("partner_electron", fYAMLConfig);
  }

  if (fConfig.ProcessDMeson()) {
    fDMesonSelection = sel::DMesonSelection(dhfe::model::kD0, fYAMLConfig);
  }
}

void AliAnalysisTaskDHFeCorr::UserCreateOutputObjects() {
  ConfigureFromYaml();
  // Additional setup for the task that needs to be done at the time of run time

  //fEventTree = std::unique_ptr<TTree>(new TTree("event", "event", 99, nullptr));
  //fElectronTree = std::unique_ptr<TTree> (new TTree("electron", "electron", 99, nullptr)); 
  //fNHFePairTree = std::unique_ptr<TTree> (new TTree("nhfe_pair", "nhfe_pair", 99, nullptr));
  //fDmesonTree = std::unique_ptr<TTree> (new TTree("dmeson", "dmeson", 99, nullptr));
  //fElectronTreeMC = std::unique_ptr<TTree>(new TTree("electron_mc", "electron_mc", 99, nullptr));
  //fDmesonTreeMC = std::unique_ptr<TTree>(new TTree("dmeson_mc", "dmeson_mc", 99, nullptr));

  if (fConfig.ProcessDMeson()) {
    fDMesonSelection.RectangularPreSelection()->GetPidHF()->SetPidResponse(
        fInputHandler->GetPIDResponse());
  }

 std::cout << "test 1" <<std::endl;
  // Remove the trigger mask from the automatic cuts
  fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny);

  fOptEvent.SetOwner(kTRUE);

  fEvent.AddToTree(*fEventTree, dhfe::configuration::kBasketSize);
  fElectron.AddToTree(*fElectronTree, dhfe::configuration::kBasketSize, IsMC());
  fDmeson.AddToTree(*fDmesonTree, dhfe::configuration::kBasketSize, IsMC());

  if (IsMC()) {
    fMCD.AddToTree(*fDmesonTreeMC, dhfe::configuration::kBasketSize);
    fMCE.AddToTree(*fElectronTreeMC, dhfe::configuration::kBasketSize);
  }

  // Event QA
  fEventCuts.AddQAplotsToList(&fOptEvent);

  // Main electron QA
  fMainESelection.AddHistogramsToOutput(fOptElectron);
  fEQAFilterBit = qa::ElectronQAHist(fMainEQAConfig, "FilterBit");
  fEQAFilterBit.AddToOutput(fOptElectron);

  fEQATrack = qa::ElectronQAHist(fMainEQAConfig, "Track");
  fEQATrack.AddToOutput(fOptElectron);

  fEQAPID = qa::ElectronQAHist(fMainEQAConfig, "PID");
  fEQAPID.AddToOutput(fOptElectron);

  // Partner electron QA
  fPartnerESelection.AddHistogramsToOutput(fOptElectron);
  fPartnerEQAFilterBit = qa::ElectronQAHist(fPartnerEQAConfig, "FilterBit");
  fPartnerEQAFilterBit.AddToOutput(fOptElectron);

  fPartnerEQATrack = qa::ElectronQAHist(fPartnerEQAConfig, "Track");
  fPartnerEQATrack.AddToOutput(fOptElectron);

  fPartnerEQAPID = qa::ElectronQAHist(fPartnerEQAConfig, "PID");
  fPartnerEQAPID.AddToOutput(fOptElectron);

  // D Meson QA
  fDMesonQAFiltering = qa::DMesonQAHist(fDMesonQAConfig, "Filtering");
  fDMesonQAFiltering.AddToOutput(fOptDMeson);

  fDMesonQAPreSelection = qa::DMesonQAHist(fDMesonQAConfig, "PreSelection");
  fDMesonQAPreSelection.AddToOutput(fOptDMeson);

  PostOutput();
}

void AliAnalysisTaskDHFeCorr::UserExec(Option_t *) {
  if (!InputEvent()) {
    PostOutput();
    return;
  }

  auto run = static_cast<unsigned int>(GetAODEvent()->GetRunNumber());
  auto event_number = static_cast<unsigned int>(fInputHandler->GetReadEntry());
  fEventId = mdl::EventId(run, std::string(CurrentFileName()), event_number);

  if (!fEventCuts.AcceptEvent(InputEvent())) {
    PostOutput();
    return;
  }

  if (fConfig.SaveEvent()) {
    fEvent = mdl::Event(fEventId, fEventCuts.GetPrimaryVertex(),
                        GetMultiSelection());
    fEventTree->Fill();
  }

  if (fConfig.ProcessElectron()) ElectronAnalysis();

  PostOutput();
}

AliAODEvent *AliAnalysisTaskDHFeCorr::GetAODEvent() const {
  return dynamic_cast<AliAODEvent *>(InputEvent());
};

void AliAnalysisTaskDHFeCorr::PostOutput() {
  PostData(1, fElectronTree.get());
  PostData(2, fDmesonTree.get());
  PostData(3, fEventTree.get());
  PostData(4, &fOptEvent);
  PostData(5, &fOptElectron);
  PostData(6, &fOptDMeson);
  PostData(7, fElectronTreeMC.get());
  PostData(8, fDmesonTreeMC.get());
}

AliMultSelection *AliAnalysisTaskDHFeCorr::GetMultiSelection() const {
  return dynamic_cast<AliMultSelection *>(
      InputEvent()->FindListObject("MultSelection"));
}

std::vector<mdl::Electron> AliAnalysisTaskDHFeCorr::SelectFilterBit(
    unsigned int filter_bit) {
  const AliAODEvent *aod_event = GetAODEvent();
  std::vector<mdl::Electron> tracks;
  tracks.reserve(static_cast<unsigned long>(aod_event->GetNumberOfTracks()));

  const TClonesArray *mc_info = GetMCInfo();
  const AliPIDResponse *pid_response = GetPIDResponse();

  for (int i(0); i < aod_event->GetNumberOfTracks(); i++) {
    auto track = dynamic_cast<AliAODTrack *>(aod_event->GetTrack(i));
    if (track->TestFilterBit(filter_bit)) {
      tracks.emplace_back(fEventId, track, i, pid_response, mc_info);
    }
  }

  tracks.shrink_to_fit();

  return tracks;
}

std::vector<mdl::Electron> AliAnalysisTaskDHFeCorr::ElectronAnalysis() {
  // Main electron selection
  auto tracks = SelectFilterBit(fMainESelection.GetFilterBit());
  fEQAFilterBit.Fill(tracks);

  auto selected_tracks = fMainESelection.FilterTracking(tracks);
  fEQATrack.Fill(selected_tracks);

  auto selected_electrons = fMainESelection.FilterPID(selected_tracks);
  fEQAPID.Fill(selected_electrons);

  // Partner electron selection
  auto partner_tracks = SelectFilterBit(fPartnerESelection.GetFilterBit());
  fPartnerEQAFilterBit.Fill(partner_tracks);

  auto selected_partner_tracks =
      fPartnerESelection.FilterTracking(partner_tracks);
  fPartnerEQATrack.Fill(selected_partner_tracks);

  auto selected_partner_electrons =
      fPartnerESelection.FilterPID(selected_partner_tracks);
  fPartnerEQATrack.Fill(selected_partner_electrons);

  // Non-HFe Reconstruction
  auto non_hfe_pairs = nhfe::FindNonHFe(
      selected_electrons, selected_partner_electrons, GetAODEvent());

  if (IsMC()) {
    auto mc = GetMCInfo();
    auto select_hfe = std::function<bool(const AliAODMCParticle *)>(
        [mc](const AliAODMCParticle *p) { return org::IsHFe(p, mc); });

    auto hfe_mc = mdl::MCParticle::MCParticles(fEventId, mc, select_hfe);
    FillTreeFromStdContainer(hfe_mc, &fMCE, *fElectronTreeMC);

    // In case of MC, it is necessary to save:
    // - All HFes in the event (even if not selected by the PID cuts)
    // - All the selected tracks (HFe or not).
    selected_electrons = AddAllHFeToTracks(selected_tracks, selected_electrons);
  }

  FillTreeFromStdContainer(selected_electrons, &fElectron, *fElectronTree);

  return selected_electrons;
}

std::vector<mdl::Electron> AliAnalysisTaskDHFeCorr::AddAllHFeToTracks(
    const std::vector<mdl::Electron> &selected_tracks,
    const std::vector<mdl::Electron> &selected_electrons) {
  std::vector<mdl::Electron> hfe_tracks;
  hfe_tracks.reserve(selected_tracks.size());

  std::function<bool(const mdl::Electron &)> is_hfe{&mdl::Electron::IsHFe};

  std::copy_if(selected_tracks.begin(), selected_tracks.end(),
               std::back_inserter(hfe_tracks), is_hfe);

  auto hfe_and_selected = std::set<mdl::Electron>(selected_electrons.begin(),
                                                  selected_electrons.end());

  hfe_and_selected.insert(hfe_tracks.begin(), hfe_tracks.end());

  std::vector<dhfe::model::Electron> all_tracks(hfe_and_selected.begin(),
                                                hfe_and_selected.end());

  return all_tracks;
}

AliPIDResponse *AliAnalysisTaskDHFeCorr::GetPIDResponse() const {
  return fInputHandler->GetPIDResponse();
}

TClonesArray *AliAnalysisTaskDHFeCorr::GetMCInfo() const {
  return dynamic_cast<TClonesArray *>(
      InputEvent()->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
}

*/