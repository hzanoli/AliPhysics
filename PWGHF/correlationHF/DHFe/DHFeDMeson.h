#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEDMESON_H_
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEDMESON_H_

#include <array>
#include <map>
#include <string>
#include <utility>

#include "AliAODRecoDecayHF.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "DHFe/DHFeEvent.h"
#include "TClonesArray.h"

namespace dhfe {
namespace model {

/* Enum with the types of D mesons supported. */
enum DSpecies_t { kD0, kDplus, kDstar };

/* Class used to store general properties of the D mesons, such as MC
 * representation, daughters, particle information. */
class DMesonData {
 public:
  /* Creates a DMesonDatabase from the provided information */
  DMesonData(std::string name, DSpecies_t species, std::string aod_list_name,
             int pdg_code, std::vector<int> daughter_pdg,
             std::vector<AliPID::EParticleType> f_daughter_ali_pid)
      : fName(std::move(name)),
        fSpecies(species),
        fAODListName(std::move(aod_list_name)),
        fPdgCode(pdg_code),
        fDaughterPDG(std::move(daughter_pdg)),
        fDaughterAliPID(std::move(f_daughter_ali_pid)) {}

 public:
  const std::string& Name() const { return fName; }
  DSpecies_t Species() const { return fSpecies; }
  const std::string& AodListName() const { return fAODListName; }
  int PdgCode() const { return fPdgCode; }
  const std::vector<int>& DaughtersPdg() const { return fDaughterPDG; }
  const std::vector<AliPID::EParticleType>& DaughtersAliPid() const {
    return fDaughterAliPID;
  }

 private:
  std::string fName;
  DSpecies_t fSpecies;
  std::string fAODListName;
  int fPdgCode;
  std::vector<int> fDaughterPDG;
  std::vector<AliPID::EParticleType> fDaughterAliPID;
};

/* Creates and stores one instance of DMesonData with the information for
 * each D meson. */
class DMesonDatabase {
 public:
  /* Given a name (D0, Dplus, Dstar), returns the corresponding DMesonData. */
  const DMesonData& FindByName(const std::string& name);

  /* Given a DSpecies_t, returns the corresponding DMesonData. */
  const DMesonData& FindSpecies(DSpecies_t species);

  /* Returns the data for the D0 meson.*/
  const DMesonData& D0() { return fD0; };

  /* Returns the data for the Dplus meson*/
  const DMesonData& Dplus() { return fDplus; };

  /*Returns the data for the D*+ meson.*/
  const DMesonData& Dstar() { return fDstar; };

 private:
  const DMesonData fD0 = DMesonData("D0", kD0, "D0toKpi", 421, {211, 321},
                                    {AliPID::kPion, AliPID::kKaon});

  const DMesonData fDplus =
      DMesonData("Dplus", kDplus, "Charm3Prong", 411, {211, 321, 211},
                 {AliPID::kPion, AliPID::kKaon, AliPID::kPion});

  const DMesonData fDstar =
      DMesonData("Dstar", kDstar, "Dstar", 413, {321, 321, 211},
                 {AliPID::kPion, AliPID::kKaon, AliPID::kPion});

  const std::vector<DMesonData> fDMesons{fD0, fDplus, fDstar};
};

/* Representation of the decay product of the D mesons */
class DecayParticle {
 public:
  /* Creates a DecayParticle by providing the information for each field. */
  DecayParticle(unsigned int id, float pt, float d0,
                AliPID::EParticleType pid_hypothesis, float tpc_response,
                float tof_response)
      : fId(id),
        fPt(pt),
        fD0(d0),
        fPIDHypothesis(pid_hypothesis),
        fTPCNSig(tpc_response),
        fTOFNSig(tof_response) {}

  /* Creates a DecayParticle by providing the mother, the id in the
   * prong track list (prong_id), the PID response task (pid_response) and
   * the particle hypothesis (pid_hypothesis). */
  DecayParticle(AliAODRecoDecayHF* mother, unsigned int prong_id,
                AliPIDResponse* pid_response,
                AliPID::EParticleType pid_hypothesis);

  /* Adds the variables of this DecayParticle to a Tree. Since it is
   * possible to have multiple instances of decay particles for the same
   * decay mother, you can use the parameter suffix (for example, 'Pion',
   * 'Pion0', 'Kaon') to add a suffix at the end of the name of the variables.*/
  void AddToTree(TTree& tree, const std::string& suffix, int buff_size);

  DecayParticle() = default;

  unsigned int Id() const { return fId; }
  float Pt() const { return fPt; }
  float D0() const { return fD0; }
  AliPID::EParticleType PIDHypothesis() const { return fPIDHypothesis; }
  float TPCNSig() const { return fTPCNSig; }
  float TOFNSig() const { return fTOFNSig; }

 private:
  unsigned int fId{0};  // ESD ID
  float fPt{-999.};
  float fD0{-999.};  // Displacement
  AliPID::EParticleType fPIDHypothesis{AliPID::kUnknown};
  float fTPCNSig{-999.};  // TPC N Sigma for the PIDHypothesis
  float fTOFNSig{-999.};  // TOF N Sigma for the PIDHypothesis
};

/* Information of a D meson in simulations. */
class DMesonMC {
 public:
  DMesonMC() = default;

  /* Construct from the candidate and MC information */
  DMesonMC(AliAODRecoDecayHF* candidate, int d_pdg,
           std::vector<int> pdg_daughter, TClonesArray* mc_info = nullptr);

  /* Adds the variables of this DMesonMC to a Tree */
  void AddToTree(TTree& tree, int buff_size);

 private:
  unsigned int fLabel{0};  // MC label of the particle
  float fPtMC{-999.};
  bool fIsD{false};         // Is it signal or background? (MC information)
  bool fIsParticle{false};  // is it a D0 or D0bar? (MC information)
  bool fIsPrompt{false};    // Does it comes from charm? (MC information)
};

class DMeson : public EventId, public DMesonMC {
 public:
  DMeson() = default;

  /* Builds a D meson from the event id, d meson id, reconstructed vertex,
   * species, pid_response and invariant mass. Avoid to use this
   * implementation. It is usually better to use one of the derived classes,
   * which fill all the fields automatically. */
  DMeson(const EventId& event_id, unsigned int d_id,
         AliAODRecoDecayHF* candidate, DSpecies_t species, float inv_mass,
         float y, int d_pdg, std::vector<int> pdg_daughter,
         TClonesArray* mc_info = nullptr);

  /* Creates branches to hold the properties of the D meson id in a tree. */
  virtual void AddToTree(TTree& tree, int buff_size, bool is_mc);

  DSpecies_t Species() const { return fSpecies; }
  unsigned int Id() const { return fID; }
  bool IsParticleCandidate() const { return fIsParticleCandidate; }
  float Pt() const { return fPt; }
  float Eta() const { return fEta; }
  float Phi() const { return fPhi; }
  float Y() const { return fY; }
  float InvMass() const { return fInvMass; }
  float ReducedChi2() const { return fReducedChi2; }
  float DecayLength() const { return fDecayLength; }
  float DecayLengthXy() const { return fDecayLengthXY; }
  float NormDecayLength() const { return fNormDecayLength; }
  float NormDecayLengthXy() const { return fNormDecayLengthXY; }
  float CosP() const { return fCosP; }
  float CosPxy() const { return fCosPXY; }
  float ImpParXy() const { return fImpParXY; }
  float Dca() const { return fDCA; }
  float Normd0MeasMinusExp() const { return fNormd0MeasMinusExp; }
  unsigned char SelectionStatusDefaultPid() const {
    return fSelectionStatusDefaultPID;
  }

 private:
  using DMesonMC::AddToTree;
  using EventId::AddToTree;

  DSpecies_t fSpecies{kD0};
  // Id of the D meson in the AOD list
  unsigned int fID{0};

  // Particle hypotheses at reconstruction
  bool fIsParticleCandidate{true};

  float fPt{-999.};
  float fEta{-999.};
  float fPhi{-999.};
  float fY{-999.};
  float fInvMass{-999.};
  float fReducedChi2{-999.};

  float fDecayLength{-999.};
  float fDecayLengthXY{-999.};

  float fNormDecayLength{-999.};
  float fNormDecayLengthXY{-999.};

  float fCosP{-999.};
  float fCosPXY{-999.};

  float fImpParXY{-999.};
  float fDCA{-999.};

  float fNormd0MeasMinusExp{-999.};
  unsigned char fSelectionStatusDefaultPID{0};

  static float ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF* candidate);
};

/* Class to store D0 candidates. They are reconstructed in the D0 -> K- pi+
 * decay channel. */
class D0 : public DMeson {
 public:
  /* Builds a D0 from an EventId, a reconstructed HF vertex and the PID response
   * task. The last variable determines if a D0 or a D0bar (its antiparticle)
   * should be constructed. */
  D0(const EventId& event_id, unsigned int d_id, AliAODRecoDecayHF* candidate,
     AliPIDResponse* pid_response, bool anti_particle = false,
     TClonesArray* mc_info = nullptr);

  /* Factory method to create D0bar ("reflections" from the D0). */
  static D0 MakeReflection(const EventId& event_id, unsigned int d_id,
                           AliAODRecoDecayHF* candidate,
                           AliPIDResponse* pid_response,
                           TClonesArray* mc_info = nullptr) {
    return {event_id, d_id, candidate, pid_response, true, mc_info};
  }

  /* Creates branches to hold the properties of the dmeson id in a tree. */
  void AddToTree(TTree& tree, int buff_size, bool is_mc) override;

  float CosTs() const { return fCosTs; };

  const DecayParticle& Pion() const { return fPion; };
  const DecayParticle& Kaon() const { return fKaon; };
  std::vector<DecayParticle> Daughters() const { return {fPion, fKaon}; };

 private:
  float fCosTs{-999.};
  DecayParticle fPion;
  DecayParticle fKaon;
};

/* Class to store D+ candidates. They are reconstructed in the D+ -> K- pi+ pi+
 * decay channel. */
class Dplus : public DMeson {
 public:
  /* Builds a D+ meson from the event id, reconstructed vertex, species,
   * pid_response. */
  Dplus(const EventId& event_id, unsigned int d_id,
        AliAODRecoDecayHF* candidate, AliPIDResponse* pid_response,
        TClonesArray* mc_info = nullptr);

  /* Creates branches to hold the properties of the dmeson id in a tree. */
  void AddToTree(TTree& tree, int buff_size, bool is_mc) override;

  float DCA1() const { return fDCA1; }
  float DCA2() const { return fDCA2; }
  float SigmaVertex() const { return fSigmaVertex; }
  const DecayParticle& PionFirst() const { return fPionFirst; }
  const DecayParticle& Kaon() const { return fKaon; }
  const DecayParticle& PionSecond() const { return fPionSecond; }

 private:
  float fDCA1{-999.};
  float fDCA2{-999.};
  float fSigmaVertex{-999.};

  DecayParticle fPionFirst;
  DecayParticle fKaon;
  DecayParticle fPionSecond;
};

}  // namespace model
}  // namespace dhfe

#endif