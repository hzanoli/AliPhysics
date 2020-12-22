#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFENONHFE_H_
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFENONHFE_H_

#include "DHFe/DHFeCorrElectron.h"
#include "DHFe/DHFeSelectionBase.h"
#include "TTree.h"

namespace dhfe {
namespace non_hfe {

/* Class to store the information of partner electrons, when making pairs to
 * identify non-HF electrons. */
class ElectronPair {
 public:
  /* Creates it using the track properties. The main_id corresponds to the id
   * of the main electron used to calculate the pair properties.*/
  enum Sign_t { kULS = -1, kLS = 1 };
  ElectronPair(unsigned int main_id, float inv_mass, float pt,
               unsigned short int crossed_rows_tpc, Sign_t sign, float red_chi2,
               float ndf)
      : fMainEId(main_id),
        fPairSign(sign),
        fInvMass(inv_mass),
        fPt(pt),
        fNCrossedRowsTPC(crossed_rows_tpc),
        fRedChi2(red_chi2),
        fNDF(ndf){};

  ElectronPair() = default;

  /* Returns true if the main and partner electron have different charge
   * signal (+- or -+). */
  bool IsULS() const { return fPairSign == kULS; };

  /* Returns true if the main and partner electron have the same charge
   * signal (-- or ++). */
  bool IsLS() const { return fPairSign == kLS; };

  /* Returns the ESD id of the main electron. */
  unsigned int MainEId() const { return fMainEId; };

  Sign_t PairSign() const { return fPairSign; };
  float InvMass() const { return fInvMass; };
  float PtPartner() const { return fPt; };
  float NCrossedRowsTPCPartner() const { return fNCrossedRowsTPC; };
  float RedChi2() const { return fRedChi2; };
  float NDF() const { return fNDF; };

  /* Add the variables of the partner to a Tree. */
  void AddToTree(TTree &tree, int buff_size);

 private:
  unsigned int fMainEId{0};
  Sign_t fPairSign{kULS};
  float fInvMass{-1};
  float fPt{-999};
  unsigned short int fNCrossedRowsTPC{0};
  float fRedChi2{-999.};
  float fNDF{-999.};
};

/* Given a list of main electrons, a vector with the partner candidates and an
 * AOD event, calculates all the Non-HFe partners for the electrons in this
 * vector. */
std::vector<ElectronPair> FindNonHFe(
    std::vector<model::Electron> &main,
    const std::vector<model::Electron> &partners, AliAODEvent *event);

/* Given a main electron, a vector with the partner candidates and an AOD
 * event, calculates all the Non-HFe partners for this electron. */
std::vector<ElectronPair> FindNonHFePartners(
    const model::Electron &main, const std::vector<model::Electron> &partners,
    AliAODEvent *event);

/* Class to store and perform the selection of non-hfe pairs.
 * The configuration of this class should be done using a yaml file,
 * following the follow model:
 *
 * non_hfe_pairs:
 *   max_inv_mass: 0.2 # The maximum invariant mass of the electron pair
 *   min_ndf: 1 # The minimum of the Number of Degrees of Freedom (NDF)
 *   max_reduced_chi2: 3. # The maximum of the Chi2/NDF fit
 */
class NonHFePairSelection {
 public:
  NonHFePairSelection() = default;
  /* Constructs the selection for the Non-HFe pairs from a configuration yaml
   * file. */
  explicit NonHFePairSelection(const PWG::Tools::AliYAMLConfiguration &yaml);

  /* Given a vector with the pairs, returns a vector combining the result
   * (true or false) of all the selections of this class. */
  std::vector<char> Select(const std::vector<ElectronPair> &pairs) const;

  /* Given a vector of non-hfe electron pairs, returns a new vector with only
   * the pairs that fulfill the selection requirements. */
  std::vector<ElectronPair> Filter(
      const std::vector<ElectronPair> &pairs) const {
    return selection::Filter(pairs, Select(pairs));
  };

  template <typename F>
  using Feature = selection::Feature<ElectronPair, F>;

  typedef selection::Selection<ElectronPair> Selection;

 private:
  float fInvMassMax{0.2};
  float fReducedChi2Max{0.5};
  float fNDFMin{1.};
};

}  // namespace non_hfe
}  // namespace dhfe

#endif