#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_SELECTION_H
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_SELECTION_H

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TH1D.h"

namespace dhfe {
namespace selection {

/* Given an object of type T and one property (or "feature") of it (of type
 * F), returns a function that selects property in the range min <=
 * property < max.
 */
template <typename T, typename F>
std::function<bool(const T &)> Range(std::function<F(const T &)> feature, F min,
                                     F max) {
  std::function<bool(const T &x)> selection = [feature, min, max](const T &x) {
    return (feature(x) >= min) && (feature(x) < max);
  };

  return selection;
}

/* Given an object of type T and one property (or "feature") of it (of type
 * F), returns a function that selects property the range property >= min.
 */
template <typename T, typename F>
std::function<bool(const T &)> Min(std::function<F(const T &)> feature, F min) {
  std::function<bool(const T &x)> selection = [feature, min](const T &x) {
    return feature(x) >= min;
  };

  return selection;
}

/* Given an object of type T and one property (or "feature") of it (of type
 * F), returns a function that selects property the range property <
 * max.
 */
template <typename T, typename F>
std::function<bool(const T &)> Max(std::function<F(const T &)> feature, F max) {
  std::function<bool(const T &x)> selection = [feature, max](const T &x) {
    return feature(x) < max;
  };

  return selection;
}

/* Given a vector with the particles and a vector with the filters (or
 * "selections") to be applied, returns a vector combining the result of all
 * the filters (combined using AND logic). */
template <typename T>
std::vector<char> Select(
    const std::vector<T> &particles,
    const std::vector<std::function<bool(const T &)>> &filters) {
  std::vector<char> combined;
  combined.reserve(particles.size());

  for (auto &&p : particles) {
    std::vector<char> selection;
    selection.reserve(filters.size());

    for (auto &&f : filters) {
      selection.push_back(f(p));
    }

    combined.push_back(std::all_of(selection.begin(), selection.end(),
                                   [](char x) { return bool(x); }));
  }

  return combined;
}

/* Given a vector with the particles and a vector with their selection status
 * (True or False, but indicated as char because of the problem with
 * std::vector<bool>), returns a vector with the particles which the
 * selection status is True. */
template <typename T>
std::vector<T> Filter(const std::vector<T> &particles,
                      const std::vector<char> &selection) {
  if (particles.size() != selection.size()) {
    throw std::runtime_error(
        "The selection and particle vectors have different sizes.");
  }

  std::vector<T> filtered;
  auto selection_it = selection.begin();
  filtered.reserve(particles.size());

  for (auto &&p : particles) {
    if (*selection_it) filtered.push_back(p);
    selection_it++;
  }

  filtered.shrink_to_fit();

  return filtered;
}

template <typename T>
using Selection = std::function<bool(const T &)>;

template <typename T, typename F>
using Feature = std::function<F(const T &)>;

template <typename T>
class Cut {
 public:
  /* Constructs a cut from a name and a selection. Use the factory
   * methods MakeCutMin, MakeCutMax and MakeCutRange to help you when
   * defining the selection. */
  Cut(std::string name, std::function<bool(const T &)> selection,
      std::string representation = "")
      : fName(std::move(name)),
        fRepresentation(std::move(representation)),
        fSelection(selection){};

  Cut() = default;

  template <typename feature_type>
  static Cut MakeCutMin(const std::string &variable_name,
                 Feature<T, feature_type> feature, feature_type min_value);

  template <typename feature_type>
  static Cut MakeCutMax(const std::string &variable_name,
                 Feature<T, feature_type> feature, feature_type max_value);

  template <typename feature_type>
  static Cut MakeCutRange(const std::string &variable_name,
                   Feature<T, feature_type> feature, feature_type min_value,
                   feature_type max_value);

  /* Apply this cut to a single object.
  bool Select(const T &p) const { return fSelection(p); };  */

  /* Apply this cyt to a vector of objects. */
  std::vector<char> Select(const std::vector<T> &vector) const;

  const std::string &Name() const { return fName; };
  const std::string &Representation() const { return fRepresentation; };

 private:
  std::string fName{"EmptyCut"};
  std::string fRepresentation{"EmptyCut"};
  std::function<bool(const T &)> fSelection{[](const T &t) { return true; }};
};

template <typename T>
std::vector<char> Cut<T>::Select(const std::vector<T> &vector) const {
  std::vector<char> selection;
  selection.reserve(vector.size());

  for (auto &&particle : vector) {
    selection.push_back(fSelection(particle));
  }

  return selection;
}

template <typename T>
template <typename feature_type>
Cut<T> Cut<T>::MakeCutMin(const std::string &variable_name,
                          Feature<T, feature_type> feature,
                          feature_type min_value) {
  std::string cut_name = variable_name + "Min";

  auto selection = Min(feature, min_value);

  std::stringstream cut_representation;
  cut_representation << variable_name << " >= " << min_value;

  return Cut(cut_name, selection, cut_representation.str());
}

template <typename T>
template <typename feature_type>
Cut<T> Cut<T>::MakeCutMax(const std::string &variable_name,
                          Feature<T, feature_type> feature,
                          feature_type max_value) {
  std::string cut_name = variable_name + "Max";

  auto selection = Max(feature, max_value);

  std::stringstream cut_representation;
  cut_representation << variable_name << " < " << max_value;

  return Cut(cut_name, selection, cut_representation.str());
}

template <typename T>
template <typename feature_type>
Cut<T> Cut<T>::MakeCutRange(const std::string &variable_name,
                            Feature<T, feature_type> feature,
                            feature_type min_value, feature_type max_value) {
  std::string cut_name = variable_name + "Range";

  auto selection = Range(feature, min_value, max_value);

  std::stringstream cut_representation;
  cut_representation << min_value << " <= " << variable_name;
  cut_representation << " < " << max_value;

  return Cut(cut_name, selection, cut_representation.str());
}

/* Class to aggregate the cuts and apply them. */
template <typename T>
class SelectionManager {
 public:
  /* Given the particle name and a vector with the cuts, creates an instance.*/
  SelectionManager(const std::string &name, const std::vector<Cut<T>> &cuts)
      : fName{name},
        fCuts{cuts},
        fCutStatus((name + "_Selection").c_str(),
                   (name + "_Selection; Cut; Counts").c_str(), fCuts.size(), 0,
                   fCuts.size()) {
    for (int i{0}; i < fCuts.size(); i++)
      fCutStatus.GetXaxis()->SetBinLabel(i + 1, fCuts[i].Name().c_str());
  };

  SelectionManager() = default;

  const std::vector<Cut<T>> &Cuts() const { return fCuts; };
  const std::string &Name() const { return fName; };

  /* Applies the selection of this SelectionManager to particles. Returns a
   * vector of the same size as particles with the status of the selection. */
  std::vector<char> Select(std::vector<T> particles);

  /* Given a TList, adds control histograms to it. The histogram monitors how
   * many particles fulfill each of the selections. */
  void AddHistogramsToOutput(TList &list);

 private:
  std::string fName;
  std::vector<Cut<T>> fCuts;
  TH1D fCutStatus;

  /* Given the selection status of the cuts, fill the number of time each cut
   * was accepted into fCutStatus. */
  void FillCutQA(const std::vector<std::vector<char>> &selection_status);

  /* Applies the cuts defined in this manager to */
  std::vector<std::vector<char>> ApplyCuts(
      const std::vector<T> &particles) const;
  std::vector<char> CombineCuts(
      const std::vector<std::vector<char>> &selection_status) const;
};

template <typename T>
std::vector<char> SelectionManager<T>::Select(std::vector<T> particles) {
  std::vector<std::vector<char>> selection_status = ApplyCuts(particles);

  FillCutQA(selection_status);

  return CombineCuts(selection_status);
}

template <typename T>
std::vector<char> SelectionManager<T>::CombineCuts(
    const std::vector<std::vector<char>> &selection_status) const {
  std::vector<char> combined;
  combined.reserve(selection_status[0].size());

  for (int i(0); i < selection_status[0].size(); i++) {
    std::vector<char> status;
    status.reserve(selection_status.size());

    for (auto &&sel : selection_status) {
      status.push_back(sel[i]);
    }

    bool combined_status = std::all_of(status.begin(), status.end(),
                                       [](char x) { return bool(x); });
    combined.push_back(combined_status);
  }
  return combined;
}

template <typename T>
std::vector<std::vector<char>> SelectionManager<T>::ApplyCuts(
    const std::vector<T> &particles) const {
  std::vector<std::vector<char>> selection_status;
  selection_status.reserve(fCuts.size());

  for (const Cut<T> &cut : fCuts) {
    selection_status.push_back(cut.Select(particles));
  }

  return selection_status;
}

template <typename T>
void SelectionManager<T>::FillCutQA(
    const std::vector<std::vector<char>> &selection_status) {
  std::vector<int> n_selected_per_cut;
  n_selected_per_cut.reserve(selection_status.size());

  for (auto &&status : selection_status) {
    int n_selected = std::count(status.begin(), status.end(), true);
    n_selected_per_cut.push_back(n_selected);
  }

  for (int i(0); i < n_selected_per_cut.size(); i++) {
    fCutStatus.Fill(i, n_selected_per_cut[i]);
  }
}
template <typename T>
void SelectionManager<T>::AddHistogramsToOutput(TList &list) {
  list.Add(&fCutStatus);
}

}  // namespace selection
}  // namespace dhfe

template <typename T>
std::ostream &operator<<(std::ostream &os, dhfe::selection::Cut<T> const &cut) {
  return os << cut.Representation();
}

template <typename T>
std::ostream &operator<<(std::ostream &os,
                         dhfe::selection::SelectionManager<T> const &manager) {
  os << "SelectionManager: " <<manager.Name();
  for (auto &&cut : manager.Cuts()) { os << std::endl << cut; }
  return os;
}
#endif