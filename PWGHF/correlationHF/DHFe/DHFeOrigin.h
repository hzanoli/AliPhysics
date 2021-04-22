#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEORIGIN_H_
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEORIGIN_H_

#include "AliAODMCParticle.h"
#include "TClonesArray.h"

namespace dhfe {
namespace origin {

enum HFOrigin_t { kCharm = 4, kBeauty = 5, kNotHF = 0 };

/* Given the array of the mc particles and the position of a particle in this
 * array (given by label), returns the origin of this particle (kCharm, kBeauty
 * or kNotHF).
 * If label refers to an electron: it is checked if it it a direct decay of
 * an HF particle and classifies the mother into */
HFOrigin_t Origin(const TClonesArray *mc, unsigned int label);

/* Given the array of the mc particles and the position of a
 * heavy-flavour particle in this array (given by label), returns kBeauty if
 * any of the mothers is a beauty meson or baryon, and  kCharm otherwise. */
HFOrigin_t HFParticleOrigin(const TClonesArray *mc, unsigned int label);

/* Given a particle pdg code and the pdg code of its mother, returns true if
 * the particle is a heavy-flavour decay electron.*/
bool IsHFe(int pdg_code, int mother_pdg);

/* Given a AliAODMCParticle returns true if the particle is a heavy-flavour
 * decay electron.*/
bool IsHFe(const AliAODMCParticle* particle, const TClonesArray *mc);

/* Given a pdg code, returns if it corresponds to a beauty meson or baryon.*/
bool IsBeauty(int pdg);

/* Given a pdg code, returns if it corresponds to a charm meson or baryon.*/
bool IsCharm(int pdg);

/* Given a pdg code, returns if it corresponds to a charm or beauty particle
 * (meson or baryon). */
bool IsHFParticle(int pdg);

}
}

#endif