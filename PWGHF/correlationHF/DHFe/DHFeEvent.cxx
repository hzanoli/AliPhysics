#include "DHFe/DHFeEvent.h"

#include <sstream>
#include <stdexcept>

#include "AliMultSelection.h"
#include "AliVVertex.h"
#include "TTree.h"

namespace dhfe {
namespace model {
void EventId::AddToTree(TTree &tree, int basked_size) {
  tree.Branch("Run", &fRun, "Run/i", basked_size);
  tree.Branch("Directory", &fDirectory, "Directory/i", basked_size);
  tree.Branch("Event", &fEvent, "Event/i", basked_size);
}

unsigned int EventId::GetFolderNumber(const std::string &file_path) {
  std::vector<std::string> folder_list;
  std::istringstream path_istring(file_path);
  std::string current_string;

  while (std::getline(path_istring, current_string, '/')) {
    folder_list.push_back(current_string);
  }
  std::string folder("0");

  if (folder_list.size() > 2)
    folder = folder_list[folder_list.size()-2];
  else
    throw std::runtime_error("Invalid file path. The task cannot handle "
        "non grid-like paths.");

  unsigned int folder_number = std::stoi(folder);
  return folder_number;
}

Event::Event(const EventId &id, const AliVVertex *vertex,
             AliMultSelection *multi_selection)
    : EventId(id) {
  fVtxZ = vertex->GetZ();

  if (multi_selection) {
    fMultV0MPercentile = multi_selection->GetMultiplicityPercentile("V0M");

    fMultiRefMult08Percentile =
        multi_selection->GetMultiplicityPercentile("RefMult08");

    fMultiSPDTrackletsPercentile =
        multi_selection->GetMultiplicityPercentile("SPDTracklets");

    if (multi_selection->GetEstimator("V0M")) {
      fMultV0M = multi_selection->GetEstimator("V0M")->GetValue();
    }

    if (multi_selection->GetEstimator("RefMult08")) {
      fMultiRefMult08 = multi_selection->GetEstimator("RefMult08")->GetValue();
    }

    if (multi_selection->GetEstimator("SPDTracklets")) {
      fMultiSPDTracklets =
          multi_selection->GetEstimator("SPDTracklets")->GetValue();
    }
  }
}
void Event::AddToTree(TTree &tree, int basked_size) {
  EventId::AddToTree(tree, basked_size);

  tree.Branch("VtxZ", &fVtxZ, "VtxZ/F", basked_size);
  tree.Branch("MultV0M", &fMultV0M, "MultV0M/F", basked_size);
  tree.Branch("MultiRefMult08", &fMultiRefMult08, "MultiRefMult08/F",
              basked_size);
  tree.Branch("MultiSPDTracklets", &fMultiSPDTracklets, "MultiSPDTracklets/F",
              basked_size);
  tree.Branch("MultV0MPercentile", &fMultV0MPercentile, "MultV0MPercentile/F",
              basked_size);
  tree.Branch("MultiRefMult08Percentile", &fMultiRefMult08Percentile,
              "MultiRefMult08Percentile/F", basked_size);
  tree.Branch("MultiSPDTrackletsPercentile", &fMultiSPDTrackletsPercentile,
              "MultiSPDTrackletsPercentile/F", basked_size);
}
}  // namespace model
}  // namespace dhfe
