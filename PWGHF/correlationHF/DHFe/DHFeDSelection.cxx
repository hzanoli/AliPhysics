#include "DHFeDSelection.h"

#include <stdexcept>

#include "AliRDHFCuts.h"
#include "TFile.h"

namespace dhfe {
namespace selection {

DMesonSelection::DMesonSelection(model::DSpecies_t species,
                                 const PWG::Tools::AliYAMLConfiguration &yaml) {
  std::string file_name;
  std::string name(model::DMesonDatabase().FindSpecies(species).Name());

  yaml.GetProperty({name, "selection", "cut_file"}, file_name, false);
  auto file_cut = std::unique_ptr<TFile>{TFile::Open(file_name.c_str())};

  std::string cut_name;
  yaml.GetProperty({name, "selection", "cut_name"}, cut_name, true);

  if (file_cut->Get(cut_name.c_str())) {
    TObject *cuts_obj = file_cut->Get(cut_name.c_str());
    auto cuts = dynamic_cast<AliRDHFCuts *>(cuts_obj);
    fRectangularCuts = std::unique_ptr<AliRDHFCuts>(cuts);
    cuts->SaveAs(("Cuts_" + name + ".root").c_str());
  } else {
    throw std::invalid_argument(
        "Check Config file. It is not possible to set"
        " the cuts for D mesons.");
  }

  yaml.GetProperty({name, "selection", "min_pt"}, fPtMin, true);
  yaml.GetProperty({name, "selection", "max_pt_pid"}, fPtMaxPID, true);
  fUsePID = this->fRectangularCuts->GetIsUsePID();
}
}  // namespace selection
}  // namespace dhfe