#include "DHFeDMeson.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "AliAODMCParticle.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"
#include "TClonesArray.h"

namespace dhfe {
namespace model {

void DMesonMC::AddToTree(TTree &tree, int buff_size) {
  tree.Branch("Label", &fLabel, "Label/i", buff_size);
  tree.Branch("PtMC", &fPtMC, "PtMC/F", buff_size);
  tree.Branch("IsD", &fIsD, "IsD/o", buff_size);
  tree.Branch("IsParticle", &fIsParticle, "IsParticle/o", buff_size);
  tree.Branch("IsPrompt", &fIsD, "IsPrompt/o", buff_size);
}
DMesonMC::DMesonMC(AliAODRecoDecayHF *candidate, int d_pdg,
                   std::vector<int> pdg_daughter, TClonesArray *mc_info) {
  if (mc_info) {
    fLabel = candidate->MatchToMC(d_pdg, mc_info, pdg_daughter.size(),
                                  &pdg_daughter[0]);
    if (fLabel < 0) return;

    const auto mc_part = dynamic_cast<AliAODMCParticle *>(mc_info->At(fLabel));
    fPtMC = mc_part->Pt();

    fIsD = true;
    fIsParticle = mc_part->PdgCode() > 0;

    if (AliVertexingHFUtils::CheckOrigin(mc_info, mc_part, false) != 5)
      fIsPrompt = true;
  }
}

DMeson::DMeson(const EventId &event_id, unsigned int d_id,
               AliAODRecoDecayHF *candidate, DSpecies_t species, float inv_mass,
               float y, int d_pdg, std::vector<int> pdg_daughter,
               TClonesArray *mc_info)
    : EventId(event_id),
      DMesonMC(candidate, d_pdg, std::move(pdg_daughter), mc_info),
      fSpecies(species),
      fID(d_id),
      fIsParticleCandidate(candidate->GetSign() >= 0),
      fPt(candidate->Pt()),
      fEta(candidate->Eta()),
      fPhi(candidate->Phi()),
      fY(y),
      fInvMass(inv_mass),
      fReducedChi2(candidate->GetReducedChi2()),
      fDecayLength(candidate->DecayLength()),
      fDecayLengthXY(candidate->DecayLengthXY()),
      fNormDecayLength(candidate->NormalizedDecayLength()),
      fNormDecayLengthXY(candidate->NormalizedDecayLengthXY()),
      fCosP(candidate->CosPointingAngle()),
      fCosPXY(candidate->CosPointingAngleXY()),
      fImpParXY(candidate->ImpParXY()),
      fDCA(candidate->GetDCA()),
      fNormd0MeasMinusExp(ComputeMaxd0MeasMinusExp(candidate)) {}

float DMeson::ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF *candidate) {
  const float b_field = candidate->GetEvent()->GetMagneticField();
  float d_d0_max = 0;
  const unsigned int n_prongs_candidate = candidate->GetNProngs();

  for (int i(0); i < n_prongs_candidate; i++) {
    double d0_diff;
    double error_d0_diff;

    candidate->Getd0MeasMinusExpProng(i, b_field, d0_diff, error_d0_diff);

    const double norm_dd0 = d0_diff / error_d0_diff;

    if ((i == 0) || (TMath::Abs(norm_dd0) > TMath::Abs(d_d0_max)))
      d_d0_max = float(norm_dd0);
  }

  return d_d0_max;
}

void DMeson::AddToTree(TTree &tree, int buff_size, bool is_mc) {
  EventId::AddToTree(tree, buff_size);
  if (is_mc) DMesonMC::AddToTree(tree, buff_size);

  // Basic information
  tree.Branch("Pt", &fPt, "Pt/F", buff_size);
  tree.Branch("Eta", &fEta, "Eta/F", buff_size);
  tree.Branch("Phi", &fPhi, "Phi/F", buff_size);
  tree.Branch("Y", &fY, "Y/F", buff_size);
  tree.Branch("InvMass", &fInvMass, "InvMass/F", buff_size);
  tree.Branch("ReducedChi2", &fReducedChi2, "ReducedChi2/F", buff_size);

  // Topological information
  tree.Branch("DecayLength", &fDecayLength, "DecayLength/F", buff_size);
  tree.Branch("DecayLengthXY", &fDecayLengthXY, "DecayLengthXY/F", buff_size);

  tree.Branch("NormDecayLength", &fNormDecayLength, "NormDecayLength/F",
              buff_size);
  tree.Branch("NormDecayLengthXY", &fNormDecayLengthXY, "NormDecayLengthXY/F",
              buff_size);

  tree.Branch("CosP", &fCosP, "CosP/F", buff_size);
  tree.Branch("CosPXY", &fCosPXY, "CosPXY/F", buff_size);

  tree.Branch("ImpParXY", &fImpParXY, "ImpParXY/F", buff_size);
  tree.Branch("DCA", &fDCA, "DCA/F", buff_size);

  tree.Branch("Normd0MeasMinusExp", &fNormd0MeasMinusExp,
              "Normd0MeasMinusExp/F", buff_size);
}

const DMesonData &DMesonDatabase::FindByName(const std::string &name) {
  auto find_by_name = [name](const DMesonData &d) { return d.Name() == name; };
  auto it = std::find_if(fDMesons.begin(), fDMesons.end(), find_by_name);

  if (it != fDMesons.end()) {
    return *it;
  }
  std::string error_message = "It is not possible to find the meson: " + name;
  throw std::invalid_argument(error_message);
}

const DMesonData &DMesonDatabase::FindSpecies(DSpecies_t species) {
  auto find_by_species = [species](const DMesonData &d) {
    return d.Species() == species;
  };
  auto it = std::find_if(fDMesons.begin(), fDMesons.end(), find_by_species);

  if (it != fDMesons.end()) {
    return *it;
  }
  throw std::invalid_argument("It is not possible to find species specified.");
}

D0::D0(const EventId &event_id, unsigned int d_id, AliAODRecoDecayHF *candidate,
       AliPIDResponse *pid_response, bool anti_particle, TClonesArray *mc_info)
    : DMeson(
          event_id, d_id, candidate, kD0,
          anti_particle
              ? dynamic_cast<AliAODRecoDecayHF2Prong *>(candidate)
                    ->InvMassD0bar()
              : dynamic_cast<AliAODRecoDecayHF2Prong *>(candidate)->InvMassD0(),
          dynamic_cast<AliAODRecoDecayHF2Prong *>(candidate)->YD0(),
          DMesonDatabase().D0().PdgCode(), DMesonDatabase().D0().DaughtersPdg(),
          mc_info) {
  const auto d0 = dynamic_cast<AliAODRecoDecayHF2Prong *>(candidate);

  if (anti_particle) {
    fCosTs = d0->CosThetaStarD0bar();
    fPion = DecayParticle(candidate, 1, pid_response, AliPID::kPion);
    fKaon = DecayParticle(candidate, 0, pid_response, AliPID::kKaon);
  } else {
    fCosTs = d0->CosThetaStarD0();
    fPion = DecayParticle(candidate, 0, pid_response, AliPID::kPion);
    fKaon = DecayParticle(candidate, 1, pid_response, AliPID::kKaon);
  }
}
void D0::AddToTree(TTree &tree, int buff_size, bool is_mc) {
  DMeson::AddToTree(tree, buff_size, is_mc);
  tree.Branch("CosTs", &fCosTs, "CosTs/F", buff_size);

  fPion.AddToTree(tree, "Pion", buff_size);
  fKaon.AddToTree(tree, "Kaon", buff_size);
}

DecayParticle::DecayParticle(AliAODRecoDecayHF *mother, unsigned int prong_id,
                             AliPIDResponse *pid_response,
                             AliPID::EParticleType pid_hypothesis) {
  auto track = dynamic_cast<AliAODTrack *>(mother->GetDaughter(prong_id));
  fId = track->GetID();
  fPt = mother->PtProng(prong_id);
  fD0 = mother->Getd0Prong(prong_id);
  fPIDHypothesis = pid_hypothesis;
  fTPCNSig = pid_response->NumberOfSigmasTPC(track, pid_hypothesis);
  fTOFNSig = pid_response->NumberOfSigmasTOF(track, pid_hypothesis);
}

void DecayParticle::AddToTree(TTree &tree, const std::string &suffix,
                              int buff_size) {
  std::string id = "Id" + suffix;
  tree.Branch(id.c_str(), &fId, (id + "/i").c_str(), buff_size);

  std::string pt = "Pt" + suffix;
  tree.Branch(pt.c_str(), &fPt, (pt + "/F").c_str(), buff_size);

  std::string d0 = "D0" + suffix;
  tree.Branch(d0.c_str(), &fD0, (d0 + "/F").c_str(), buff_size);

  std::string tpc = "TPCNSig" + suffix;
  tree.Branch(tpc.c_str(), &fTPCNSig, (tpc + "/F").c_str(), buff_size);

  std::string tof = "TOFNSig" + suffix;
  tree.Branch(tof.c_str(), &fTPCNSig, (tof + "/F").c_str(), buff_size);
}

Dplus::Dplus(const EventId &event_id, unsigned int d_id,
             AliAODRecoDecayHF *candidate, AliPIDResponse *pid_response,
             TClonesArray *mc_info)
    : DMeson(event_id, d_id, candidate, kDplus,
             dynamic_cast<AliAODRecoDecayHF3Prong *>(candidate)->InvMassDplus(),
             dynamic_cast<AliAODRecoDecayHF3Prong *>(candidate)->YDplus(),
             DMesonDatabase().Dplus().PdgCode(),
             DMesonDatabase().Dplus().DaughtersPdg(), mc_info) {
  fDCA1 = candidate->GetDCA(1);
  fDCA2 = candidate->GetDCA(2);

  const auto dplus = dynamic_cast<AliAODRecoDecayHF3Prong *>(candidate);
  fSigmaVertex = dplus->GetSigmaVert();

  fPionFirst = DecayParticle(candidate, 0, pid_response, AliPID::kPion);
  fKaon = DecayParticle(candidate, 1, pid_response, AliPID::kKaon);
  fPionSecond = DecayParticle(candidate, 2, pid_response, AliPID::kPion);
}

void Dplus::AddToTree(TTree &tree, int buff_size, bool is_mc) {
  DMeson::AddToTree(tree, buff_size, is_mc);
  fPionFirst.AddToTree(tree, "Pion0", buff_size);
  fKaon.AddToTree(tree, "Kaon", buff_size);
  fPionSecond.AddToTree(tree, "Pion1", buff_size);

  tree.Branch("DCA1", &fDCA1, "DCA1/F", buff_size);
  tree.Branch("DCA2", &fDCA2, "DCA2/F", buff_size);
  tree.Branch("SigmaVertex", &fSigmaVertex, "SigmaVertex/F", buff_size);
}
}  // namespace model
}  // namespace dhfe