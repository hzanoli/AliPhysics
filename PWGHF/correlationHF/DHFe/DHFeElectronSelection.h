#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEELECTRONSELECTION_H_
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEELECTRONSELECTION_H_
#include "AliAODTrack.h"
#include "AliYAMLConfiguration.h"
#include "DHFeCorrElectron.h"
#include "DHFeSelectionBase.h"

namespace dhfe {
namespace selection {

enum ITSPixel_t {
  kFirst = 0,            // at least one on first
  kSecond = 1,           // at least one on second
  kBoth = 2,             // hit on layer 1 and 2
  kNone = 3,             // no hits
  kAny = 4,              // have at least one hit on layer 1 or 2
  kExclusiveSecond = 5,  // only on layer 2
  kExclusiveFirst = 6    // only on layer 1
};

class ElectronSelectionYamlMappings {
 public:
  const std::map<std::string, AliAODTrack::AODTrkFilterBits_t> &AODFilterBits()
      const {
    return fgkAODFilterBitMap;
  }
  const std::map<std::string, ITSPixel_t> &ITSPixel() const {
    return fgkITSPixelMap;
  }
  const std::string &GetITSPixelString(ITSPixel_t pixel_req);

 private:
  const std::map<std::string, AliAODTrack::AODTrkFilterBits_t>
      fgkAODFilterBitMap = {
          {"kTrkTPCOnly", AliAODTrack::kTrkTPCOnly},
          {"kTrkITSsa", AliAODTrack::kTrkITSsa},
          {"kTrkITSConstrained", AliAODTrack::kTrkITSConstrained},
          {"kTrkElectronsPID", AliAODTrack::kTrkElectronsPID},
          {"kTrkGlobalNoDCA", AliAODTrack::kTrkGlobalNoDCA},
          {"kTrkGlobal", AliAODTrack::kTrkGlobal},
          {"kTrkGlobalSDD", AliAODTrack::kTrkGlobalSDD},
          {"kTrkTPCOnlyConstrained", AliAODTrack::kTrkTPCOnlyConstrained}};
  const std::map<std::string, ITSPixel_t> fgkITSPixelMap = {
      {"kFirst", kFirst},
      {"kSecond", kSecond},
      {"kBoth", kBoth},
      {"kNone", kNone},
      {"kAny", kAny},
      {"kExclusiveSecond", kExclusiveSecond},
      {"kExclusiveFirst", kExclusiveFirst}};
};

/* Class to store and perform the selection of electrons.
  The configuration of this class should be done using a yaml file,
  following the follow model:

particle: # Name of the particle, such as main_electron, electron, etc
  track:
    pt_range: [minimum, maximum]
    eta_range: [minimum, maximum]
    filter_bit: kTrkGlobalNoDCA # Check possible in
        # ElectronSelectionYamlMappings::fgkAODFilterBitMap
    tpc_crossed_rows_min: 0
    tpc_cls_dedx_min: 0
    its_hits_min: 2
    pixel_req: kAny # check possibilities in ITSPixel_t
    dca_z_max: 5.0
    dca_xy_max: 5.0
  PID:
    tpc_nsigma_range: [minimum, maximum]
    require_tof: false
    tof_nsigma_range: [minimum, maximum]
*/
class ElectronSelection {
 public:
  ElectronSelection() = default;

  /* Constructs the selection for a particle (main_electron or
   * partner_electron) from a configuration yaml file */
  ElectronSelection(const std::string &particle,
                    const PWG::Tools::AliYAMLConfiguration &config);

  /* Template definition of classes to use as Feature for electrons*/
  template <typename F>
  using Feature = selection::Feature<model::Electron, F>;

  AliAODTrack::AODTrkFilterBits_t GetFilterBit() const { return fFilterBit; }

  /* Given a vector of electrons, returns a new vector with only the electrons
   * that fulfill the tracking requirements.*/
  std::vector<model::Electron> FilterTracking(
      const std::vector<model::Electron> &electrons) {
    return selection::Filter(electrons, SelectTracking(electrons));
  };

  /* Given a vector of electrons, returns a new vector with only the electrons
   * that fulfill the PID requirements.*/
  std::vector<model::Electron> FilterPID(
      const std::vector<model::Electron> &electrons) {
    return selection::Filter(electrons, SelectPID(electrons));
  };

  /* Given an electron e and a requirement on the ITS detector (one of the cases
   * defined ITSPixel_t), returns true if this e fulfills this requirement. */
  static bool FulfilPixelSelection(const model::Electron &e,
                                   ITSPixel_t requirement);

  /* Given a TList, adds control histograms to it. The histogram monitors how
   * many particles fulfill each of the selections. */
  void AddHistogramsToOutput(TList &list);

  const std::string& ParticleName() const { return fParticleName; };

  const SelectionManager<dhfe::model::Electron>
      &TrackSelectionManager() const {
    return fTrackCuts;
  };

  const SelectionManager<dhfe::model::Electron>
      &PIDSelectionManager() const {
    return fPIDCuts;
  };

 private:
  /* Given a vector of electrons, performs the tracking selection (Pt, Eta,
   * NCrossedRowsTPC, TPCClsDeDx, ITSCls, ITSPixel, DCAxy and DCAz). Returns a
   * vector with the selection status (True or False) for each electron. */
  std::vector<char> SelectTracking(
      const std::vector<model::Electron> &electrons);

  /* Given a vector of electrons, performs the tracking selection (Only TPC and
   * TOF for electron hypotheses). Returns a vector with the selection status
   * (True or False) for each electron. */
  std::vector<char> SelectPID(const std::vector<model::Electron> &electrons);

  /* Reads the track selection from the yaml configuration. */
  void ReadTrackCuts(const std::string &particle,
                     const PWG::Tools::AliYAMLConfiguration &config);

  /* Reads the PID selection from the yaml configuration. */
  void ReadPIDCuts(const std::string &particle,
                   const PWG::Tools::AliYAMLConfiguration &config);

  /* Creates a selection manager from the track selections defined in this
   * class.*/
  SelectionManager<dhfe::model::Electron> CreateTrackCuts();

  /* Creates a selection manager from the PID selections defined in this
   * class.*/
  SelectionManager<dhfe::model::Electron> CreatePIDCuts();

  std::string fParticleName;
  SelectionManager<dhfe::model::Electron> fTrackCuts;
  SelectionManager<dhfe::model::Electron> fPIDCuts;

  // Track selection
  Float_t fPtMin{Float_t(-999.)};
  Float_t fPtMax{Float_t(999.)};
  Float_t fEtaMin{Float_t(-999.)};
  Float_t fEtaMax{Float_t(999.)};

  AliAODTrack::AODTrkFilterBits_t fFilterBit{AliAODTrack::kTrkGlobalNoDCA};

  UShort_t fNCrossedRowsTPCMin{0};
  UShort_t fTPCClsDeDxMin{0};
  UChar_t fITSClsMin{0};

  ITSPixel_t fITSPixel{kAny};

  Float_t fDCAzMax{Float_t(999.)};
  Float_t fDCAxyMax{Float_t(999.)};

  // PID selection
  Bool_t fRequireTOF{kFALSE};
  Float_t fTOFNSigmaMin{Float_t(-3.)};
  Float_t fTOFNSigmaMax{Float_t(3.)};

  Float_t fTPCNSigmaElectronMin{Float_t(-3.)};
  Float_t fTPCNSigmaElectronMax{Float_t(3.)};
};

}  // namespace selection
}  // namespace dhfe

std::ostream &operator<<(std::ostream &os,
                         dhfe::selection::ElectronSelection const
                             &e_selection);

#endif
