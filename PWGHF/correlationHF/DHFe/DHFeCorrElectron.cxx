#include "DHFeCorrElectron.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVertexingHFUtils.h"
#include "DHFeOrigin.h"
#include "TClonesArray.h"

namespace sel = dhfe::selection;
namespace config = dhfe::configuration;

namespace dhfe {
namespace model {

Track::Track(const AliAODTrack *track, unsigned int id_aod) {
  fId = track->GetID();
  fAodId = id_aod;
  fCharge = track->Charge();
  fPt = track->Pt();
  fP = track->P();
  fEta = track->Eta();
  fPhi = track->Phi();

  fNCrossedRowsTPC = track->GetTPCNCrossedRows();
  fNClsTPCDeDx = track->GetTPCsignalN();
  fNITSCls = track->GetITSNcls();

  fITSHitFirstLayer = track->HasPointOnITSLayer(0);
  fITSHitSecondLayer = track->HasPointOnITSLayer(1);

  Double_t d0z0[2] = {-999., -999.};
  Double_t cov[3] = {-999., -999., -999.};

  const auto event = track->GetAODEvent();
  const AliVVertex *primaryVertex = event->GetPrimaryVertex();

  if (AliAODTrack(*track).PropagateToDCA(
          primaryVertex, event->GetMagneticField(), 20., d0z0, cov)) {
    fDCAxy = d0z0[0];
    fDCAz = d0z0[1];
  }
}

void Track::AddToTree(TTree &tree, int buff_size) {
  tree.Branch("Id", &fId, "Id/i", buff_size);
  tree.Branch("Pt", &fPt, "Pt/F", buff_size);
  tree.Branch("Eta", &fEta, "Eta/F", buff_size);
  tree.Branch("Phi", &fPhi, "Phi/F", buff_size);
  tree.Branch("Charge", &fCharge, "Charge/B", buff_size);
  tree.Branch("P", &fP, "P/F", buff_size);
  tree.Branch("NCrossedRowsTPC", &fNCrossedRowsTPC, "fNCrossedRowsTPC/s",
              buff_size);
  tree.Branch("NClsTPCDeDx", &fNClsTPCDeDx, "NClsTPCDeDx/s", buff_size);
  tree.Branch("NITSCls", &fNITSCls, "NITSCls/b", buff_size);
  tree.Branch("ITSHitFirstLayer", &fITSHitFirstLayer, "ITSHitFirsMakeBinstLayer/O",
              buff_size);
  tree.Branch("ITSHitSecondLayer", &fITSHitSecondLayer, "ITSHitSecondLayer/O",
              buff_size);
  tree.Branch("DCAxy", &fDCAxy, "DCAxy/F", buff_size);
  tree.Branch("DCAz", &fDCAz, "DCAz/F", buff_size);
}

ParticleIdentification::ParticleIdentification(const AliAODTrack *track,
                                               const AliPIDResponse *pid) {
  fTPCNSigmaElectron = pid->NumberOfSigmasTPC(track, AliPID::kElectron);
  fTPCNSigmaProton = pid->NumberOfSigmasTPC(track, AliPID::kProton);
  fTPCNSigmaKaon = pid->NumberOfSigmasTPC(track, AliPID::kKaon);
  fTPCNSigmaPion = pid->NumberOfSigmasTPC(track, AliPID::kPion);

  fTOFNSigmaElectron = pid->NumberOfSigmasTOF(track, AliPID::kElectron);
  fTOFNSigmaProton = pid->NumberOfSigmasTOF(track, AliPID::kProton);
  fTOFNSigmaKaon = pid->NumberOfSigmasTOF(track, AliPID::kKaon);
  fTOFNSigmaPion = pid->NumberOfSigmasTOF(track, AliPID::kPion);
}

void ParticleIdentification::AddToTree(TTree &tree, int buff_size) {
  tree.Branch("TPCNSigmaElectron", &fTPCNSigmaElectron, "TPCNSigmaElectron/F",
              buff_size);
  tree.Branch("TPCNSigmaProton", &fTPCNSigmaProton, "TPCNSigmaProton/F",
              buff_size);
  tree.Branch("TPCNSigmaKaon", &fTPCNSigmaKaon, "TPCNSigmaKaon/F", buff_size);
  tree.Branch("TPCNSigmaPion", &fTPCNSigmaPion, "TPCNSigmaPion/F", buff_size);

  tree.Branch("TOFNSigmaElectron", &fTOFNSigmaElectron, "TOFNSigmaElectron/F",
              buff_size);
  tree.Branch("TOFNSigmaProton", &fTOFNSigmaProton, "TOFNSigmaProton/F",
              buff_size);
  tree.Branch("TOFNSigmaKaon", &fTOFNSigmaKaon, "TOFNSigmaKaon/F", buff_size);
  tree.Branch("TOFNSigmaPion", &fTOFNSigmaPion, "TOFNSigmaPion/F", buff_size);
}

RecoElectronMC::RecoElectronMC(int label, const TClonesArray *mc_info) {
  if (mc_info) {
    fIsFilled = true;
    fLabel = label;
    auto part = dynamic_cast<AliAODMCParticle *>(mc_info->At(abs(label)));

    if (part) {
      fPtMC = part->Pt();
      fEtaMC = part->Eta();
      fPhiMC = part->Phi();
      fPDGCode = part->PdgCode();
      fOrigin = dhfe::origin::Origin(mc_info, label);

      if (part->GetMother() > 0) {
        const auto mother =
            dynamic_cast<AliAODMCParticle *>(mc_info->At(part->GetMother()));
        fFirstMotherPDG = mother->GetPdgCode();
        fFirstMotherPt = mother->Pt();

        if (mother->GetMother() > 0) {
          const auto g_mother = dynamic_cast<AliAODMCParticle *>(
              mc_info->At(mother->GetMother()));
          fSecondMotherPDG = g_mother->GetPdgCode();
          fSecondMotherPt = g_mother->Pt();
        }
      }
    }
  }
}

void RecoElectronMC::AddToTree(TTree &tree, int buff_size) {
  tree.Branch("PtMC", &fPtMC, "PtMC/F", buff_size);
  tree.Branch("Label", &fLabel, "Label/I", buff_size);
  tree.Branch("Origin", &fOrigin, "Origin/b", buff_size);
  tree.Branch("PDGCode", &fPDGCode, "PDGCode/I", buff_size);
  tree.Branch("FirstMotherPDG", &fFirstMotherPDG, "FirstMotherPDG/I",
              buff_size);
  tree.Branch("FirstMotherPt", &fFirstMotherPt, "FirstMotherPt/F", buff_size);
  tree.Branch("SecondMotherPDG", &fSecondMotherPDG, "SecondMotherPDG/I",
              buff_size);
  tree.Branch("SecondMotherPt", &fSecondMotherPt, "SecondMotherPt/F",
              buff_size);
}

bool RecoElectronMC::IsHFe() const {
  return dhfe::origin::IsHFe(fPDGCode, fFirstMotherPDG);
}

void Electron::AddToTree(TTree &tree, int buff_size, bool is_mc) {
  EventId::AddToTree(tree, buff_size);
  Track::AddToTree(tree, buff_size);
  ParticleIdentification::AddToTree(tree, buff_size);
  if (is_mc) RecoElectronMC::AddToTree(tree, buff_size);
}

}  // namespace model
}  // namespace dhfe
