#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFECORRQA_H_
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFECORRQA_H_

#include <utility>
#include <vector>

#include "AliYAMLConfiguration.h"
#include "DHFe/DHFeCorrElectron.h"
#include "DHFe/DHFeDMeson.h"
#include "DHFe/DHFeNonHFe.h"
#include "DHFe/DHFeUtils.h"

namespace dhfe {
namespace qa {
/* Returns evenly spaced bins in the interval [start,end] with n_bins. */
std::vector<Float_t> MakeBins(Float_t start, Float_t end, Int_t n_bins);

/* Configures the Electron QA binning */
class ElectronQAConfig {
 public:
  ElectronQAConfig() = default;
  /* Given the name of the particle (usually MainE or PartnerE) and a file with
   * the yaml configuration, creates the QA configuration. */
  ElectronQAConfig(const std::string &name,
                   const PWG::Tools::AliYAMLConfiguration &yaml_config);

  const std::string &ParticleName() const { return fParticleName; }
  const std::vector<Float_t> &PtBins() const { return fPtBins; }
  const std::vector<Float_t> &PBins() const { return fPBins; }
  Int_t NBinsEta() const { return fNBinsEta; }
  Int_t NBinsPhi() const { return fNBinsPhi; }
  Int_t NBinsTpcCls() const { return fNBinsTPCCls; }
  Int_t NBinsNsigma() const { return fNBinsNsigma; }

 private:
  std::string fParticleName{"MainE"};
  std::vector<Float_t> fPtBins{0, 0.5, 1.0, 1.5, 2.0, 3., 4., 6., 10.};
  std::vector<Float_t> fPBins{0, 0.5, 1.0, 1.5, 2.0, 3., 4., 6., 10.};
  Int_t fNBinsEta{20};
  Int_t fNBinsPhi{40};
  Int_t fNBinsTPCCls{80};
  Int_t fNBinsNsigma{100};
};

/* Base class for QA histograms */
class QAHist {
 public:
  QAHist() = default;
  QAHist(std::string particle, std::string stage)
      : fParticle(std::move(particle)), fStage(std::move(stage)) {}

 public:
  std::string GetName(const std::string &name) {
    return (fParticle + "_" + name + "_" + fStage);
  };

 private:
  std::string fParticle;
  std::string fStage;
};

/* Class to hold and create the QA plots for the electrons.*/
class ElectronQAHist : public QAHist {
 public:
  ElectronQAHist() = default;
  /* Constructs from a configuration object and stage (for example filterbit,
   * tracking, pid).*/
  ElectronQAHist(const ElectronQAConfig &qa_config, const std::string &stage);

  /* Fill the histograms of the QA histograms with the values from electron. */
  void Fill(const model::Electron &electron);

  /* Fill the histograms of this QA object with the values for each electron. */
  void Fill(const std::vector<model::Electron> &electrons);

  void AddToOutput(TList& list);

  const TH3F &PtEtaPhi() const { return fPtEtaPhi; }
  const TH1F &TPCNCrossedRows() const { return fTPCNCrossedRows; }
  const TH1F &TPCNClsDeDx() const { return fTPCNClsDeDx; }
  const TH1F &ItsCls() const { return fITSCls; }
  const TH1F &DCAz() const { return fDCAz; }
  const TH1F &DCAxy() const { return fDCAxy; }
  const TH2F &TPCNsigmaPt() const { return fTPCNsigmaPt; }
  const TH2F &TPCNsigmaP() const { return fTPCNsigmaP; }
  const TH2F &TOFNsigmaPt() const { return fTOFNsigmaPt; }
  const TH2F &TOFNsigmaP() const { return fTOFNsigmaP; }
  const TH2F &TPCTOFNSigma() const { return fTPCTOFNSigma; }

 private:
  TH3F fPtEtaPhi;
  TH1F fTPCNCrossedRows;
  TH1F fTPCNClsDeDx;
  TH1F fITSFirstLayerHit;
  TH1F fITSSecondLayerHit;
  TH1F fITSCls;
  TH1F fDCAz;
  TH1F fDCAxy;
  TH2F fTPCNsigmaPt;
  TH2F fTPCNsigmaP;
  TH2F fTOFNsigmaPt;
  TH2F fTOFNsigmaP;
  TH2F fTPCTOFNSigma;
};

/* Class to hold and create the configuration for the QA plots for the D
 * mesons.*/
class DMesonQAConfig {
 public:
  DMesonQAConfig() = default;
  /* Constructs from a stage (for example filter, Preselection) and a yaml file
   * object.*/
  DMesonQAConfig(const std::string &meson,
                 const PWG::Tools::AliYAMLConfiguration &yaml_config);

  const std::string &ParticleName() const { return fMeson; };
  const std::string &DMeson() const { return fMeson; };
  const std::vector<Float_t> &PtBins() const { return fPtBins; }
  Float_t InvMassMin() const { return fInvMassMin; }
  Float_t InvMassMax() const { return fInvMassMax; }
  Int_t NBinsInvMass() const { return fNBinsInvMass; }
  Int_t NBinsPhi() const { return fNBinsPhi; }
  Int_t NBinsEta() const { return fNBinsEta; }

 private:
  std::string fMeson{"D0"};
  std::vector<Float_t> fPtBins{0, 1., 2., 4., 6., 8., 12., 16., 24., 36., 50.};
  Float_t fInvMassMin{1.3};
  Float_t fInvMassMax{2.5};
  Int_t fNBinsInvMass{100};

  Int_t fNBinsPhi{40};
  Int_t fNBinsEta{16};
};

/* Class to hold and create the QA plots for D mesons.*/
class DMesonQAHist : public QAHist {
 public:
  DMesonQAHist() = default;
  DMesonQAHist(const DMesonQAConfig &config, std::string stage);

  void AddToOutput(TList& list);
 private:
  TH3F fPtEtaPhi;
  TH2F fPtInvMass;
};

}  // namespace qa
}  // namespace dhfe

#endif