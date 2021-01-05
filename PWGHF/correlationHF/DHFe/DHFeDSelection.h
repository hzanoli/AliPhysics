#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEDSELECTION_H_
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEDSELECTION_H_

#include "AliYAMLConfiguration.h"
#include "DHFeDMeson.h"
#include "DHFeSelectionBase.h"
#include "AliRDHFCutsD0toKpi.h"

namespace dhfe {
namespace selection {

/* Class to store and perform the selection of D mesons.
  The configuration of this class should be done using a yaml file,
  following the follow model:

D0: # Name of the particle (D0, Dplus, Dstar)
  selection:
    cut_file: "PreselectionDHFe.root"
    cut_name: "D0toKpiCuts"
    min_pt: 2.0
    max_pt_pid: 5.0
*/
class DMesonSelection {
 public:
  DMesonSelection() = default;
  DMesonSelection(model::DSpecies_t species,
                  const PWG::Tools::AliYAMLConfiguration& yaml);

  const std::unique_ptr<AliRDHFCuts>& RectangularPreSelection() {
    return fRectangularCuts;
  };

 private:
  std::unique_ptr<AliRDHFCuts> fRectangularCuts{new AliRDHFCutsD0toKpi()};
  float fPtMin{0.};
  float fPtMaxPID{10.};
  bool fUsePID{true};
};

}  // namespace selection
}  // namespace dhfe
#endif