#include "DHFeCorrQA.h"

#include <utility>

#include "TMath.h"

namespace mdl = dhfe::model;
namespace qa = dhfe::qa;
namespace yaml = dhfe::yaml;

std::vector<Float_t> qa::MakeBins(Float_t start, Float_t end, Int_t n_bins) {
  const auto step = (end - start) / n_bins;
  std::vector<Float_t> bins;
  bins.reserve(n_bins + 1);

  for (Int_t i(0); i < (n_bins + 1); i++) bins.push_back(start + i * step);

  return bins;
}

qa::ElectronQAConfig::ElectronQAConfig(
    const std::string &name,
    const PWG::Tools::AliYAMLConfiguration &yaml_config) {
  std::string base = "qa";

  yaml_config.GetProperty({base, name, "name"}, fParticleName, true);
  yaml_config.GetProperty({base, name, "pt_bins"}, fPtBins, true);
  yaml_config.GetProperty({base, name, "p_bins"}, fPBins, true);
  yaml_config.GetProperty({base, name, "n_eta_bins"}, fNBinsEta, true);
  yaml_config.GetProperty({base, name, "n_phi_bins"}, fNBinsPhi, true);
  yaml_config.GetProperty({base, name, "n_bins_TPC_cls"}, fNBinsTPCCls, true);
  yaml_config.GetProperty({base, name, "n_bins_n_sigma"}, fNBinsNsigma, true);
}

qa::ElectronQAHist::ElectronQAHist(const ElectronQAConfig &qa_config,
                                   const std::string &stage)
    : QAHist(qa_config.ParticleName(), stage) {
  std::vector<Float_t> phi_bins =
      MakeBins(0., 2. * TMath::Pi(), qa_config.NBinsPhi());
  std::vector<Float_t> eta_bins = MakeBins(-0.8, 0.8, qa_config.NBinsEta());

  fPtEtaPhi = TH3F(GetName("PtEtaPhi").c_str(),
                   (GetName("PtEtaPhi") + ";p_{T};#eta;#varphi").c_str(),
                   qa_config.PtBins().size() - 1, &qa_config.PtBins()[0],
                   qa_config.NBinsEta(), &eta_bins[0], qa_config.NBinsPhi(),
                   &phi_bins[0]);

  fTPCNCrossedRows = TH1F(GetName("TPCNCls").c_str(),
                          (GetName("TPCNCls") + ";N cluster TPC").c_str(),
                          qa_config.NBinsTpcCls(), 0, 160);

  fTPCNClsDeDx = TH1F(GetName("TPCNClsDeDx").c_str(),
                      (GetName("TPCNClsDeDx") + ";N cluster TPC dE/dx").c_str(),
                      qa_config.NBinsTpcCls(), 0, 160);

  fITSCls = TH1F(GetName("ITSCls").c_str(),
                 (GetName("ITSCls") + ";N ITS Cls").c_str(), 7, 0, 7);

  fITSFirstLayerHit =
      TH1F(GetName("ITSFirstLayerHit").c_str(),
           (GetName("ITSFirstLayerHit") + "; Hit?; Count").c_str(), 2, 0, 2);
  fITSFirstLayerHit.GetXaxis()->SetBinLabel(1, "False");
  fITSFirstLayerHit.GetXaxis()->SetBinLabel(2, "True");

  fITSSecondLayerHit =
      TH1F(GetName("ITSSecondLayerHit").c_str(),
           (GetName("ITSSecondLayerHit") + "; Hit?; Counts").c_str(), 2, -0.001,
           1.999);
  fITSSecondLayerHit.GetXaxis()->SetBinLabel(1, "False");
  fITSSecondLayerHit.GetXaxis()->SetBinLabel(2, "True");

  fDCAz = TH1F(GetName("DCAz").c_str(), (GetName("DCAz") + ";DCA z").c_str(),
               25, 0, 5);

  fDCAxy = TH1F(GetName("DCAxy").c_str(),
                (GetName("DCAxy") + "; DCA xy").c_str(), 25, 0, 5);

  std::vector<Float_t> n_sigma_bins =
      MakeBins(-10, 10, qa_config.NBinsNsigma());

  fTPCNsigmaPt = TH2F(GetName("TPCNsigmaPt").c_str(),
                      (GetName("TPCNsigmaPt") + ";p_{T};N#sigma TPC").c_str(),
                      qa_config.PtBins().size() - 1, &qa_config.PtBins()[0],
                      qa_config.NBinsNsigma(), &n_sigma_bins[0]);

  fTPCNsigmaP = TH2F(GetName("TPCNsigmaPt").c_str(),
                     (GetName("TPCNsigmaPt") + ";p;N#sigma TPC").c_str(),
                     qa_config.PBins().size() - 1, &qa_config.PBins()[0],
                     qa_config.NBinsNsigma(), &n_sigma_bins[0]);

  fTOFNsigmaPt = TH2F(GetName("TOFNsigmaPt").c_str(),
                      (GetName("TOFNsigmaPt") + ";p_{T};N#sigma TOF").c_str(),
                      qa_config.PBins().size() - 1, &qa_config.PBins()[0],
                      qa_config.NBinsNsigma(), &n_sigma_bins[0]);

  fTOFNsigmaP = TH2F(GetName("TOFNsigmaP").c_str(),
                     (GetName("TOFNsigmaP") + ";p;N#sigma TOF").c_str(),
                     qa_config.PBins().size() - 1, &qa_config.PBins()[0],
                     qa_config.NBinsNsigma(), &n_sigma_bins[0]);

  fTPCTOFNSigma =
      TH2F(GetName("TPCTOFNSigma").c_str(),
           (GetName("TPCTOFNSigma") + ";N#sigma TPC;N#sigma TOF").c_str(),
           qa_config.NBinsNsigma(), &n_sigma_bins[0], qa_config.NBinsNsigma(),
           &n_sigma_bins[0]);
}

void qa::ElectronQAHist::Fill(const mdl::Electron &electron) {
  fPtEtaPhi.Fill(electron.Pt(), electron.Eta(), electron.Phi());
  fTPCNCrossedRows.Fill(electron.NCrossedRowsTPC());
  fTPCNClsDeDx.Fill(electron.NClsTPCDeDx());
  fITSCls.Fill(electron.NITSCls());
  fITSFirstLayerHit.Fill(electron.ITSHitFirstLayer());
  fITSSecondLayerHit.Fill(electron.ITSHitSecondLayer());
  fDCAz.Fill(electron.DCAz());
  fDCAxy.Fill(electron.DCAxy());
  fTPCNsigmaPt.Fill(electron.Pt(), electron.TPCNSigmaElectron());
  fTPCNsigmaP.Fill(electron.P(), electron.TPCNSigmaElectron());
  fTOFNsigmaPt.Fill(electron.Pt(), electron.TOFNSigmaElectron());
  fTOFNsigmaP.Fill(electron.P(), electron.TOFNSigmaElectron());
  fTPCTOFNSigma.Fill(electron.TPCNSigmaElectron(),
                     electron.TOFNSigmaElectron());
}

void qa::ElectronQAHist::Fill(const std::vector<mdl::Electron> &electrons) {
  for (auto &&e : electrons) Fill(e);
}
void dhfe::qa::ElectronQAHist::AddToOutput(TList &list) {
  list.Add(&fPtEtaPhi);
  list.Add(&fTPCNCrossedRows);
  list.Add(&fTPCNClsDeDx);
  list.Add(&fITSCls);
  list.Add(&fITSFirstLayerHit);
  list.Add(&fITSSecondLayerHit);
  list.Add(&fDCAz);
  list.Add(&fDCAxy);
  list.Add(&fTPCNsigmaPt);
  list.Add(&fTPCNsigmaP);
  list.Add(&fTOFNsigmaPt);
  list.Add(&fTOFNsigmaP);
  list.Add(&fTPCTOFNSigma);
}

qa::DMesonQAConfig::DMesonQAConfig(
    const std::string &meson,
    const PWG::Tools::AliYAMLConfiguration &yaml_config) {
  yaml::GetPropertyRange<float>(yaml_config,
                                {meson, "selection", "inv_mass_range"},
                                fInvMassMin, fInvMassMax);
  yaml_config.GetProperty({meson, "selection", "pt_bins"}, fPtBins, true);
  yaml_config.GetProperty({meson, "selection", "n_bins_inv_mass"},
                          fNBinsInvMass, true);
  yaml_config.GetProperty({meson, "selection", "n_eta_bins"}, fNBinsEta, true);
  yaml_config.GetProperty({meson, "selection", "n_phi_bins"}, fNBinsPhi, true);
}

qa::DMesonQAHist::DMesonQAHist(const DMesonQAConfig &config, std::string stage)
    : QAHist(config.ParticleName(), std::move(stage)) {
  std::vector<Float_t> phi_bins =
      MakeBins(0.0, 2 * TMath::Pi(), config.NBinsPhi());
  std::vector<Float_t> eta_bins = MakeBins(-0.8, 0.8, config.NBinsEta());

  fPtEtaPhi =
      TH3F(GetName("PtEtaPhi").c_str(),
           (GetName("PtEtaPhi") + ";p_{T};#eta;#varphi").c_str(),
           config.PtBins().size() - 1, &config.PtBins()[0], config.NBinsEta(),
           &eta_bins[0], config.NBinsPhi(), &phi_bins[0]);

  std::vector<Float_t> inv_mass_bins =
      MakeBins(config.InvMassMin(), config.InvMassMax(), config.NBinsInvMass());

  fPtInvMass = TH2F(GetName("PtInvMass").c_str(),
                    (GetName("PtInvMass") + ";p_{T};mass").c_str(),
                    config.PtBins().size() - 1, &config.PtBins()[0],
                    config.NBinsInvMass(), &inv_mass_bins[0]);
}
void dhfe::qa::DMesonQAHist::AddToOutput(TList &list) {
  list.Add(&fPtEtaPhi);
  list.Add(&fPtInvMass);
}
