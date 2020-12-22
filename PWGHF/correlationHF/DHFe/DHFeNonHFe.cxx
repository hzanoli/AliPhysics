#include "DHFeNonHFe.h"

#include <string>

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliKFParticle.h"
#include "DHFeUtils.h"

namespace sel = dhfe::selection;

namespace dhfe {
namespace non_hfe {
void ElectronPair::AddToTree(TTree &tree, int buff_size) {
  tree.Branch("ElectronId", &fMainEId, "ElectronId/i", buff_size);
  tree.Branch("Pt", &fPt, "Eta/F", buff_size);
  tree.Branch("InvMass", &fInvMass, "Eta/F", buff_size);
  tree.Branch("NCrossedRowsTPC", &fNCrossedRowsTPC, "fNCrossedRowsTPC/s",
              buff_size);
}

std::vector<ElectronPair> FindNonHFe(
    std::vector<model::Electron> &main,
    const std::vector<model::Electron> &partners, AliAODEvent *event) {

  std::vector<ElectronPair> non_hfe_pairs;
  non_hfe_pairs.reserve(partners.size());

  for (auto && e: main){
    auto pairs = FindNonHFePartners(e, partners,event);
    non_hfe_pairs.insert(non_hfe_pairs.end(), pairs.begin(), pairs.end());
  }

  non_hfe_pairs.shrink_to_fit();
  return non_hfe_pairs;
}

std::vector<ElectronPair> FindNonHFePartners(
    const model::Electron &main, const std::vector<model::Electron> &partners,
    AliAODEvent *event) {
  std::vector<ElectronPair> non_hfe_partners;
  non_hfe_partners.reserve(partners.size());

  AliKFParticle::SetField(event->GetMagneticField());

  const auto track_main = main.AODTrack(event);
  const int pdg_main = -11 * yaml::sgn(main.Charge());

  auto vtrack_main = static_cast<AliVTrack *>(track_main);
  const auto kf_main = AliKFParticle(*vtrack_main, pdg_main);

  for (const auto &partner : partners) {
    const auto track_partner = partner.AODTrack(event);

    if (TMath::Abs(track_partner->GetID()) == TMath::Abs(track_main->GetID()))
      continue;

    const int pdg_partner = -11 * yaml::sgn(partner.Charge());

    const auto vtrack_partner = static_cast<AliVTrack *>(track_partner);
    const auto kf_partner = AliKFParticle(*vtrack_partner, pdg_partner);

    AliKFParticle photon(kf_main, kf_partner);

    ElectronPair::Sign_t sign(ElectronPair::kULS);
    if (main.Charge() * partner.Charge() > 0) sign = ElectronPair::kLS;

    ElectronPair pair(
        main.Id(), photon.GetMass(), partner.Pt(), partner.NCrossedRowsTPC(),
        sign, float(photon.GetChi2() / photon.GetNDF()), photon.GetNDF());

    non_hfe_partners.push_back(pair);
  }

  non_hfe_partners.shrink_to_fit();
  return non_hfe_partners;
}

NonHFePairSelection::NonHFePairSelection(
    const PWG::Tools::AliYAMLConfiguration &yaml) {
  typedef std::vector<std::string> vector_string;

  auto inv_mass = vector_string({"non_hfe_pairs", "max_inv_mass"});
  yaml.GetProperty(inv_mass, fInvMassMax, true);

  auto red_chi2 = vector_string({"non_hfe_pairs", "max_reduced_chi2"});
  yaml.GetProperty(red_chi2, fReducedChi2Max, true);

  auto ndf = vector_string({"non_hfe_pairs", "min_ndf"});
  yaml.GetProperty(ndf, fNDFMin, true);
}
std::vector<char> NonHFePairSelection::Select(
    const std::vector<ElectronPair> &pairs) const {

  auto inv_mass_func = Feature<float>(&ElectronPair::InvMass);
  auto inv_mass_max = sel::Max(inv_mass_func, fInvMassMax);

  auto red_chi2_func = Feature<float>(&ElectronPair::RedChi2);
  auto red_chi2_max = sel::Max(red_chi2_func, fReducedChi2Max);

  auto ndf_func = Feature<float>(&ElectronPair::NDF);
  auto ndf_min = sel::Min(ndf_func, fNDFMin);

  std::vector<Selection> selections = {ndf_min, red_chi2_max, inv_mass_max};

  return sel::Select(pairs, selections);
}
}  // namespace non_hfe
}  // namespace dhfe
