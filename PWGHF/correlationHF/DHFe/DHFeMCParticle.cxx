#include "DHFe/DHFeMCParticle.h"

#include "DHFe/DHFeOrigin.h"

namespace dhfe {
namespace model {

MCParticle::MCParticle(const EventId& id, const TClonesArray* mc,
                       unsigned int label)
    : EventId(id) {
  const auto particle = dynamic_cast<AliAODMCParticle*>(mc->At(label));
  fLabel = label;
  fOrigin = dhfe::origin::Origin(mc, label);

  fE = particle->E();
  fPt = particle->Pt();
  fEta = particle->Eta();
  fPhi = particle->Phi();
  fXv = particle->Xv();
  fYv = particle->Yv();
  fZv = particle->Zv();
  fTv = particle->Tv();
  fCharge = particle->Charge();
  fPDGCode = particle->PdgCode();
}

void MCParticle::AddToTree(TTree &tree, int buff_size) {
  EventId::AddToTree(tree, buff_size);
  tree.Branch("Label", &fLabel, "Label/i", buff_size);
  tree.Branch("E", &fE, "Float/F", buff_size);
  tree.Branch("Pt", &fPt, "Pt/F", buff_size);
  tree.Branch("Eta", &fEta, "Eta/F", buff_size);
  tree.Branch("Phi", &fPhi, "Phi/F", buff_size);
  tree.Branch("Xv", &fXv, "Xv/F", buff_size);
  tree.Branch("Yv", &fYv, "Yv/F", buff_size);
  tree.Branch("Zv", &fZv, "Zv/F", buff_size);
  tree.Branch("Tv", &fTv, "Tv/F", buff_size);
  tree.Branch("Charge", &fCharge, "Charge/F", buff_size);
  tree.Branch("PDGCode", &fPDGCode, "PDGCode/F", buff_size);
  tree.Branch("Origin", &fOrigin, "Origin/F", buff_size);
}

std::vector<MCParticle> MCParticle::MCParticles(
    const EventId &id, const TClonesArray *mc,
    const std::function<bool(const AliAODMCParticle *)>& condition_to_fill) {
  std::vector<MCParticle> particles;
  particles.reserve(mc->GetEntriesFast());

  for (int i(0); i < mc->GetEntriesFast(); i++) {
    auto part = dynamic_cast<AliAODMCParticle*>(mc->At(i));

    if (!part)
      continue;

    if (condition_to_fill(part))
      particles.emplace_back(id, mc, i);
  }

  return particles;
}

}  // namespace model
}  // namespace dhfe