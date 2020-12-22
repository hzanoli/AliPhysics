#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEMCPARTICLE_H_
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEMCPARTICLE_H_

#include <functional>

#include "AliAODMCParticle.h"
#include "DHFeEvent.h"

namespace dhfe {
namespace model {

/* Stores MC information for generated particles */
class MCParticle : public EventId {
 public:
  MCParticle() = default;

  /* Creates a MC particle. */
  MCParticle(const EventId& id, const TClonesArray* mc,
             unsigned int label);

  static std::vector<MCParticle> MCParticles(
      const EventId& id, const TClonesArray* mc,
      const std::function<bool(const AliAODMCParticle*)>& condition_to_fill);

  /* Creates branches to hold the properties of the MC particle id in a tree. */
  void AddToTree(TTree& tree, int buff_size) override;

  unsigned int Label() const { return fLabel; }
  float E() const { return fE; }
  float Pt() const { return fPt; }
  float Eta() const { return fEta; }
  float Phi() const { return fPhi; }
  float Xv() const { return fXv; }
  float Yv() const { return fYv; }
  float Zv() const { return fZv; }
  float Tv() const { return fTv; }
  char Charge() const { return fCharge; }
  int PdgCode() const { return fPDGCode; }
  unsigned short Origin() const { return fOrigin; }

 private:
  unsigned int fLabel{0};

  float fE{-999.};
  float fPt{-999.};
  float fEta{-999.};
  float fPhi{-999.};

  float fXv{-999.};
  float fYv{-999.};
  float fZv{-999.};
  float fTv{-999.};

  char fCharge{99};
  int fPDGCode{-999};
  unsigned short fOrigin{99};
};
}  // namespace model
}  // namespace dhfe
#endif
