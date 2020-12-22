#include "DHFeElectronSelection.h"

#include <iostream>
#include <stdexcept>

#include "DHFeUtils.h"

namespace sel = dhfe::selection;
namespace mdl = dhfe::model;

namespace dhfe {
namespace selection {
ElectronSelection::ElectronSelection(
    const std::string &particle, const PWG::Tools::AliYAMLConfiguration &config)
    : fParticleName(particle) {
  std::cout << "Configuring track and PID for: " << particle << std::endl;

  ReadTrackCuts(particle, config);

  ReadPIDCuts(particle, config);

  fTrackCuts = CreateTrackCuts();
  fPIDCuts = CreatePIDCuts();
}
void ElectronSelection::ReadPIDCuts(
    const std::string &particle,
    const PWG::Tools::AliYAMLConfiguration &config) {
  yaml::GetPropertyRange(config, {particle, "PID", "tpc_nsigma_range"},
                         fTPCNSigmaElectronMin, fTPCNSigmaElectronMax);
  config.GetProperty({particle, "PID", "require_tof"}, fRequireTOF, true);
  if (fRequireTOF) {
    yaml::GetPropertyRange(config, {particle, "PID", "tof_nsigma_range"},
                           fTOFNSigmaMin, fTOFNSigmaMax);
  }
}

void ElectronSelection::ReadTrackCuts(
    const std::string &particle,
    const PWG::Tools::AliYAMLConfiguration &config) {
  yaml::GetPropertyRange(config, {particle, "track", "pt_range"}, fPtMin,
                         fPtMax);
  yaml::GetPropertyRange(config, {particle, "track", "eta_range"}, fEtaMin,
                         fEtaMax);

  yaml::GetPropertyMap(config, {particle, "track", "filter_bit"},
                       ElectronSelectionYamlMappings().AODFilterBits(),
                       fFilterBit);

  config.GetProperty({particle, "track", "tpc_crossed_rows_min"},
                     fNCrossedRowsTPCMin, true);
  config.GetProperty({particle, "track", "tpc_cls_dedx_min"}, fTPCClsDeDxMin,
                     true);
  config.GetProperty({particle, "track", "its_hits_min"}, fITSClsMin, true);
  yaml::GetPropertyMap(config, {particle, "track", "pixel_req"},
                       ElectronSelectionYamlMappings().ITSPixel(), fITSPixel);

  config.GetProperty({particle, "track", "dca_z_max"}, fDCAzMax, true);
  config.GetProperty({particle, "track", "dca_xy_max"}, fDCAxyMax, true);
}

sel::SelectionManager<mdl::Electron> ElectronSelection::CreateTrackCuts() {
  typedef sel::Cut<mdl::Electron> ECut;
  Feature<Float_t> pt(&model::Electron::Pt);
  ECut pt_range = ECut::MakeCutRange("Pt", pt, fPtMin, fPtMax);

  Feature<Float_t> eta(&model::Electron::Eta);
  ECut eta_range = ECut::MakeCutRange("Eta", eta, fEtaMin, fEtaMax);

  Feature<UShort_t> crossed_rows_tpc(&model::Electron::NCrossedRowsTPC);
  ECut crossed_rows_tpc_min =
      ECut::MakeCutMin("CrossedRowTPC", crossed_rows_tpc, fNCrossedRowsTPCMin);

  Feature<UShort_t> cls_tpc_dedx(&model::Electron::NClsTPCDeDx);
  ECut cls_tpc_dedx_min =
      ECut::MakeCutMin("ClsTPCDeDx", cls_tpc_dedx, fTPCClsDeDxMin);

  Feature<UChar_t> its_clusters(&model::Electron::NITSCls);
  ECut its_cls_min = ECut::MakeCutMin("ITSCls", its_clusters, fITSClsMin);

  ITSPixel_t pixel = fITSPixel;
  Selection<model::Electron> its_pixel_sel = [pixel](const model::Electron &e) {
    return FulfilPixelSelection(e, pixel);
  };
  ECut its_pixel(
      "ITSPixel", its_pixel_sel,
      "ITSPixel == " +
          ElectronSelectionYamlMappings().GetITSPixelString(fITSPixel));

  Feature<Float_t> dca_xy(&model::Electron::DCAxy);
  ECut dca_xy_max = ECut::MakeCutMax("DCAXY", dca_xy, fDCAxyMax);

  Feature<Float_t> dca_z(&model::Electron::DCAz);
  ECut dca_z_max = ECut::MakeCutMax("DCAZ", dca_z, fDCAzMax);

  std::vector<ECut> track_cuts = {
      pt_range, eta_range, crossed_rows_tpc_min, cls_tpc_dedx_min,
      its_pixel, dca_xy_max, dca_z_max}; //its_cls_min

  return {fParticleName + "_Track", track_cuts};
}

sel::SelectionManager<mdl::Electron> ElectronSelection::CreatePIDCuts() {
  typedef sel::Cut<mdl::Electron> ECut;
  std::vector<ECut> pid_cuts;

  Feature<Float_t> tpc_nsigma_e(&model::Electron::TPCNSigmaElectron);
  ECut tpc = ECut::MakeCutRange("TPCNSigmaE", tpc_nsigma_e,
                                fTPCNSigmaElectronMin, fTPCNSigmaElectronMax);
  pid_cuts.push_back(tpc);

  if (fRequireTOF) {
    Feature<Float_t> tof_nsigma_e(&model::Electron::TOFNSigmaElectron);
    ECut tof = ECut::MakeCutRange("TOFNSigmaE", tof_nsigma_e, fTOFNSigmaMin,
                                  fTOFNSigmaMax);

    pid_cuts.push_back(tof);
  }
  return {fParticleName + "_PID", pid_cuts};
}

std::vector<char> ElectronSelection::SelectTracking(
    const std::vector<model::Electron> &electrons) {
  return fTrackCuts.Select(electrons);
}
std::vector<char> ElectronSelection::SelectPID(
    const std::vector<model::Electron> &electrons) {
  return fPIDCuts.Select(electrons);
}

bool ElectronSelection::FulfilPixelSelection(const model::Electron &e,
                                             const ITSPixel_t requirement) {
  const bool first = e.ITSHitFirstLayer();
  const bool second = e.ITSHitSecondLayer();

  switch (requirement) {
    case kAny: {
      return (first || second);
    }
    case kBoth: {
      return (first && second);
    }
    case kFirst: {
      return (first);
    }
    case kSecond: {
      return (second);
    }
    case kNone: {
      return (!first && !second);
    }
    case kExclusiveFirst: {
      return (first && !second);
    }
    case kExclusiveSecond: {
      return (!first && second);
    }
  }

  throw std::runtime_error("ElectronSelection::FulfilPixelSelection: the pixed"
      " requirement is not recognized.");
}
void ElectronSelection::AddHistogramsToOutput(TList &list) {
  fTrackCuts.AddHistogramsToOutput(list);
  fPIDCuts.AddHistogramsToOutput(list);
}
const std::string &ElectronSelectionYamlMappings::GetITSPixelString(
    ITSPixel_t pixel_req) {
  for (auto &&item : fgkITSPixelMap) {
    if (item.second == pixel_req) return item.first;
  }

  throw std::runtime_error("The requested ITS pixel is not defined.");
}
}  // namespace selection
}  // namespace dhfe

std::ostream &operator<<(
    std::ostream &os, dhfe::selection::ElectronSelection const &e_selection) {
  os << "ElectronSelection: " << e_selection.ParticleName() << std::endl;
  os << e_selection.TrackSelectionManager() << std::endl;
  os << e_selection.PIDSelectionManager();

  return os;
}