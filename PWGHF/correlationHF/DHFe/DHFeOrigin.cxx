#include "DHFeOrigin.h"

namespace dhfe {
namespace origin {

HFOrigin_t Origin(const TClonesArray *mc, unsigned int label) {
  auto particle = dynamic_cast<AliAODMCParticle *>(mc->At(label));
  if (!particle) return kNotHF;

  auto mother = dynamic_cast<AliAODMCParticle *>(mc->At(particle->GetMother()));
  if (!mother) return kNotHF;

  if (abs(particle->GetPdgCode()) == 11) {
    if (!IsHFe(particle->GetPdgCode(), mother->GetPdgCode())) return kNotHF;
    return HFParticleOrigin(mc, label);
  }

  if (IsHFParticle(particle->GetPdgCode())) {
    return HFParticleOrigin(mc, label);
  }

  return kNotHF;
}

HFOrigin_t HFParticleOrigin(const TClonesArray *mc, unsigned int label) {
  const auto particle = dynamic_cast<AliAODMCParticle *>(mc->At(label));
  auto mother = dynamic_cast<AliAODMCParticle *>(mc->At(particle->GetMother()));

  while (mother) {
    if (IsBeauty(mother->PdgCode())) return kBeauty;
    mother = dynamic_cast<AliAODMCParticle *>(mc->At(mother->GetMother()));
  }
  return kCharm;
}

bool IsHFe(int pdg_code, int mother_pdg) {
  if (abs(pdg_code) != 11) return false;
  return IsHFParticle(mother_pdg);
}


bool IsHFe(const AliAODMCParticle* particle, const TClonesArray *mc) {
  if (!particle) return false;

  auto mother = dynamic_cast<AliAODMCParticle *>(mc->At(particle->GetMother()));
  if (!mother) return false;

  int pdg_code = particle->GetPdgCode();
  int mother_pdg =mother->GetPdgCode();

  if (abs(pdg_code) != 11) return false;

  return IsHFParticle(mother_pdg);
}

bool IsBeauty(int pdg) {
  pdg = abs(pdg);
  return ((pdg == 5) || (pdg > 5000 && pdg < 6000) || (pdg > 500 && pdg < 600));
}

bool IsCharm(int pdg) {
  pdg = abs(pdg);
  return ((pdg == 4) || (pdg > 4000 && pdg < 5000) || (pdg > 400 && pdg < 500));
}

bool IsHFParticle(int pdg) {
  pdg = abs(pdg);
  return (IsBeauty(pdg) || IsCharm(pdg));
}

}  // namespace origin
}  // namespace dhfe
