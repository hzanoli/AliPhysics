#ifndef AliAnalysisTaskDHFeCorr_H
#define AliAnalysisTaskDHFeCorr_H

#include <memory>
#include <string>
#include <vector>

#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPID.h"
#include "AliRDHFCuts.h"
#include "AliYAMLConfiguration.h"


#include "DHFeCorrElectron.h"
#include "DHFeCorrQA.h"
#include "DHFeDMeson.h"
#include "DHFeDSelection.h"
#include "DHFeElectronSelection.h"
#include "DHFeMCParticle.h"
#include "DHFeNonHFe.h"


#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"


namespace dhfe {
namespace config {
/* Handles the configuration of the main task. It is constructed using the YAML
configuration file that should have the following format:

task:
  mc_mode: false # Sets the task to work with simulated data
  qa_histograms: true # Keeps the qa histograms
  only_efficiency: false # Ignore background for D mesons to save space (MC only)
  process_electron: true # whether to perform the electron analysis
  process_dmeson: true # whether to perform the dmeson analysis
  save_event: true # whether to save the event information
*/    
class DHFeTaskConfig {   
 public:
  DHFeTaskConfig() = default;

  explicit DHFeTaskConfig(const PWG::Tools::AliYAMLConfiguration &config);

  bool IsMC() const { return fIsMC; }
  bool SaveHistograms() const { return fSaveHistograms; };
  bool CalculateOnlyEfficiency() const { return fIsEffMode; };
  bool ProcessElectron() const { return fProcessElectron; };
  bool ProcessDMeson() const { return fProcessDMeson; };
  bool SaveEvent() const { return fSaveEvent; };

 private:
  bool fIsMC{false};
  bool fSaveHistograms{false};
  bool fIsEffMode{false};
  bool fProcessElectron{true};
  bool fProcessDMeson{true};
  bool fSaveEvent{true};
};
}  // namespace config
}  // namespace dhfe

class AliAnalysisTaskDHFeCorr : public AliAnalysisTaskSE {
 public:
  // D-meson species
  typedef enum {
    kNotSelected = 0,
    kParticle = 1,
    kAntiParticle = 2,
    kAmbiguous = 3  // Selected for both
  } SelectionStatus_t;

  AliAnalysisTaskDHFeCorr();

  explicit AliAnalysisTaskDHFeCorr(
      const char *name, unsigned int trigger = AliVEvent::kINT7,
      const char *config =
          "$ALICE_PHYSICS/PWGHF/correlationHF/macros/default_config_d_hfe.yaml",
      const char *default_config =
          "$ALICE_PHYSICS/PWGHF/correlationHF/macros/default_config_d_hfe.yaml");

  // Connects to the AliAnalysisManager and adds the task to it
  void ConnectToAnalysisManager();

  // Given the name of the task, creates and connects the input and output
  void CreateAndConnectIO(const std::string &name);

  ~AliAnalysisTaskDHFeCorr();

  void UserCreateOutputObjects(); 
  /*
  void UserExec(Option_t *option);

  //bool IsMC() const { return fConfig.IsMC(); }
    */

 private:
   void DefineIO();

  // Configuration of the task
  std::string fDefaultConfig{"$ALICE_PHYSICS/PWGHF/correlationHF/macros/default_config_d_hfe.yaml"}; 
  std::string fUserConfig{"$ALICE_PHYSICS/PWGHF/correlationHF/macros/default_config_d_hfe.yaml"};

  // Output variables that will be saved to the ROOT file
  TList fOptEvent; //!
  TList fOptElectron; //!
  TList fOptDMeson; //!

  std::unique_ptr<TTree> fEventTree; //!
  std::unique_ptr<TTree> fElectronTree; //!
  std::unique_ptr<TTree> fNHFePairTree; //!
  std::unique_ptr<TTree> fDmesonTree; //!
  std::unique_ptr<TTree> fElectronTreeMC; //!
  std::unique_ptr<TTree> fDmesonTreeMC; //!

  // Posts the output of the task to the output containers.
  void PostOutput();

  // Selection classes for the physics analysis
  AliEventCuts fEventCuts;

  dhfe::config::DHFeTaskConfig fConfig;  //!

  /*
  PWG::Tools::AliYAMLConfiguration fYAMLConfig; 
  dhfe::selection::ElectronSelection fMainESelection;     //!
  dhfe::selection::ElectronSelection fPartnerESelection;  //!
  dhfe::selection::DMesonSelection fDMesonSelection;      //!

  // Main and partner electron QA histogram configuration
  dhfe::qa::ElectronQAConfig fMainEQAConfig;     //!
  dhfe::qa::ElectronQAConfig fPartnerEQAConfig;  //!

  // Non-HFe tagging selection
  dhfe::non_hfe::NonHFePairSelection fPhotonSelection;  //!

  // D meson QA histogram configuration
  dhfe::qa::DMesonQAConfig fDMesonQAConfig;  //!

  dhfe::model::DMesonDatabase fDMesonDatabase;  //!



  // Objects to be saved in the tree
  dhfe::model::EventId fEventId;    //!
  dhfe::model::Event fEvent;        //!
  dhfe::model::Electron fElectron;  //!
  dhfe::model::DMeson fDmeson;      //!
  dhfe::model::MCParticle fMCE;     //!
  dhfe::model::MCParticle fMCD;     //!

  

  // Main electron QA after filter bit, track and PID
  dhfe::qa::ElectronQAHist fEQAFilterBit;  //!
  dhfe::qa::ElectronQAHist fEQATrack;      //!
  dhfe::qa::ElectronQAHist fEQAPID;        //!

  // Partner electron QA after filter bit, track and PID
  dhfe::qa::ElectronQAHist fPartnerEQAFilterBit;  //!
  dhfe::qa::ElectronQAHist fPartnerEQATrack;      //!
  dhfe::qa::ElectronQAHist fPartnerEQAPID;        //!

  // D meson QA after the filtering and preselection cuts
  dhfe::qa::DMesonQAHist fDMesonQAFiltering;     //!
  dhfe::qa::DMesonQAHist fDMesonQAPreSelection;  //!

  AliAODEvent *GetAODEvent() const;

  // Uses the Yaml file configuration to load the parameters for the task.
  void ConfigureFromYaml();

  std::vector<dhfe::model::Electron> ElectronAnalysis();

  // Given a vector with the tracks that fulfill the track cuts
  // (selected_tracks) and another one that fulfills both track and PID
  // selection (selected_electrons), returns a vector with the content of
  // selected_electrons plus the HFe tracks in selected_tracks.
  static std::vector<dhfe::model::Electron> AddAllHFeToTracks(
      const std::vector<dhfe::model::Electron> &selected_tracks,
      const std::vector<dhfe::model::Electron> &selected_electrons);

  template <typename T>
  void FillTreeFromStdContainer(std::vector<T> &items, T *item_to_fill,
                                TTree &tree);

  // Returns the multiplicity selection task.
  AliMultSelection *GetMultiSelection() const;

  // Returns the array with the generated MC particles. 
  TClonesArray *GetMCInfo() const;

  // Returns the PID response task. 
  AliPIDResponse *GetPIDResponse() const;

  // Given an AOD filter bit, returns the tracks in this event that fulfill
  // this filter bit selection.
  std::vector<dhfe::model::Electron> SelectFilterBit(unsigned int filter_bit);

  */
  ClassDef(AliAnalysisTaskDHFeCorr, 6);
  
};

/*
template <typename T>
void AliAnalysisTaskDHFeCorr::FillTreeFromStdContainer(std::vector<T> &items,
                                                       T *item_to_fill,
                                                       TTree &tree) {
  for (const auto &item : items) {
    *item_to_fill = item;
    tree.Fill();
  }
}
*/


#endif
