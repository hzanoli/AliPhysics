#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFECORRELECTRON_H_
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFECORRELECTRON_H_

#include <vector>

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliYAMLConfiguration.h"
#include "DHFe/DHFeConfig.h"
#include "DHFe/DHFeEvent.h"
#include "DHFe/DHFeOrigin.h"
#include "DHFe/DHFeSelectionBase.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"

namespace dhfe {
namespace model {

/* Track: stores the properties that depend only on the track. */
class Track {
 public:
  /* Creates using an AliAODTrack. */
  explicit Track(const AliAODTrack *track, unsigned int id_aod);
  Track() = default;

  /* Creates branches to hold the properties of the track in a tree. */
  void AddToTree(TTree &tree, int buff_size);

  AliAODTrack *AODTrack(const AliAODEvent *event) const {
    return dynamic_cast<AliAODTrack *>(event->GetTrack(fAodId));
  }

  unsigned int Id() const { return fId; }
  char Charge() const { return fCharge; }
  float Pt() const { return fPt; }
  float P() const { return fP; }
  float Eta() const { return fEta; }
  float Phi() const { return fPhi; }

  unsigned short int NCrossedRowsTPC() const { return fNCrossedRowsTPC; }
  unsigned short int NClsTPCDeDx() const { return fNClsTPCDeDx; }
  unsigned char NITSCls() const { return fNITSCls; }
  bool ITSHitFirstLayer() const { return fITSHitFirstLayer; }
  bool ITSHitSecondLayer() const { return fITSHitSecondLayer; }
  float DCAxy() const { return fDCAxy; }
  float DCAz() const { return fDCAz; }

 private:
  unsigned int fId{0};     // ESD ID (unique for each track)
  unsigned int fAodId{0};  // AOD Id, to get it from the event
  char fCharge{0};
  float fPt{-999.};
  float fP{-999.};
  float fEta{-999.};
  float fPhi{-999.};

  unsigned short int fNCrossedRowsTPC{0};
  unsigned short int fNClsTPCDeDx{0};
  unsigned char fNITSCls{0};
  bool fITSHitFirstLayer{false};
  bool fITSHitSecondLayer{false};
  float fDCAxy{999.};
  float fDCAz{999.};
};

/* PID information. Stores the PID response for different detectors. */
class ParticleIdentification {
 public:
  ParticleIdentification() = default;
  /* Creates from a track and the PID response object */
  ParticleIdentification(const AliAODTrack *track, const AliPIDResponse *pid);

  /* Creates branches to hold the PID information in a tree. */
  void AddToTree(TTree &tree, int buff_size);

  float TPCNSigmaElectron() const { return fTPCNSigmaElectron; }
  float TPCNSigmaProton() const { return fTPCNSigmaProton; }
  float TPCNSigmaKaon() const { return fTPCNSigmaKaon; }
  float TPCNSigmaPion() const { return fTPCNSigmaPion; }
  float TOFNSigmaElectron() const { return fTOFNSigmaElectron; }
  float TOFNSigmaProton() const { return fTOFNSigmaProton; }
  float TOFNSigmaKaon() const { return fTOFNSigmaKaon; }
  float TOFSigmaPion() const { return fTOFNSigmaPion; }

 private:
  float fTPCNSigmaElectron{-999.};
  float fTPCNSigmaProton{-999.};
  float fTPCNSigmaKaon{-999.};
  float fTPCNSigmaPion{-999.};

  float fTOFNSigmaElectron{-999.};
  float fTOFNSigmaProton{-999.};
  float fTOFNSigmaKaon{-999.};
  float fTOFNSigmaPion{-999.};
};

/* Reconstructed electron MC information */
class RecoElectronMC {
 public:
  /* Creates from a MC label and an MC array */
  RecoElectronMC(int label, const TClonesArray *mc_info);
  RecoElectronMC() = default;

  /* Creates branches to hold the properties of the track in a tree. */
  void AddToTree(TTree &tree, int buff_size);

  /* Returns true if the MC information has been filled, false if not */
  bool IsFilled() const { return fIsFilled; }

  unsigned int Label() const { return fLabel; }
  float PtMC() const { return fPtMC; }
  float PhiMC() const { return fPhiMC; }
  float EtaMC() const { return fEtaMC; }
  unsigned short int Origin() const { return fOrigin; }
  int PDGCode() const { return fPDGCode; }
  int FirstMotherPDG() const { return fFirstMotherPDG; }
  float FirstMotherPt() const { return fFirstMotherPt; }
  int SecondMotherPDG() const { return fSecondMotherPDG; }
  float SecondMotherPt() const { return fSecondMotherPt; }

  /* Returns true if the electron comes the direct decay of a charm or beauty
   * hadrons. */
  bool IsHFe() const;

 private:
  bool fIsFilled{false};

  unsigned int fLabel{0};

  float fPtMC{-999.};
  float fPhiMC{-999.};
  float fEtaMC{-999.};

  dhfe::origin::HFOrigin_t fOrigin{dhfe::origin::HFOrigin_t::kNotHF};
  int fPDGCode{0};
  int fFirstMotherPDG{0};
  float fFirstMotherPt{-999.};
  int fSecondMotherPDG{0};
  float fSecondMotherPt{-999.};
};

/* Class to store an electron candidate. */
class Electron : public EventId,
                 public Track,
                 public ParticleIdentification,
                 public RecoElectronMC {
 public:
  Electron() = default;

  /* Creates an electron candidate from a track and event-related information.
   * If in simulations, you can pass a TClonesArray with the MC particles to
   * fill the MC information. */
  Electron(const EventId &event_id, const AliAODTrack *track,
           unsigned int id_aod, const AliPIDResponse *pid_response,
           const TClonesArray *mc_info = nullptr)
      : EventId(event_id),
        Track(track, id_aod),
        ParticleIdentification(track, pid_response),
        RecoElectronMC(abs(track->GetLabel()), mc_info){};

  bool operator==(const Electron &rhs) const { return Id() == rhs.Id(); }
  bool operator!=(const Electron &rhs) const { return !(rhs == *this); };

  bool operator<(const Electron &rhs) const { return Id() < rhs.Id(); }
  bool operator>(const Electron &rhs) const { return rhs < *this; }
  bool operator<=(const Electron &rhs) const { return !(rhs < *this); }
  bool operator>=(const Electron &rhs) const { return !(*this < rhs); }

  /* Creates a branch to hold the properties of reconstructed electrons in a
   * tree. */
  void AddToTree(TTree &tree, int buff_size, bool is_mc);

 private:
  using EventId::AddToTree;
  using ParticleIdentification::AddToTree;
  using Track::AddToTree;
};
}  // namespace model
}  // namespace dhfe

#endif
