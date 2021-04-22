#ifndef ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEUTILS_H_
#define ALIPHYSICS_PWGHF_CORRELATIONHF_DHFE_DHFEUTILS_H_

#include <exception>
#include <map>
#include <string>
#include <vector>

#include "AliYAMLConfiguration.h"

namespace dhfe {
namespace yaml {

/* Given a vector with the path to the property in the yaml file and the
 * exception raised, constructs an error message to be shown. */
std::string BuildErrorMessagePath(const std::vector<std::string> &property_path,
                                  const std::exception &exp);

/* Reads a property with minimum (property_min) and maximum (property_max) range
 * in the path property_path from a yaml file represented by yaml. */
template <typename T>
void GetPropertyRange(const PWG::Tools::AliYAMLConfiguration &yaml,
                      const std::vector<std::string> &property_path,
                      T &property_min, T &property_max) {
  std::vector<T> range;
  yaml.GetProperty(property_path, range, true);

  try {
    property_min = range[0];
    property_max = range[1];
  } catch (std::exception &exp) {
    throw std::invalid_argument(BuildErrorMessagePath(property_path, exp));
  }
}

/* Reads a property which can take values given by property_map in the path
 * property_path from a yaml file represented by yaml. */
template <typename T>
void GetPropertyMap(const PWG::Tools::AliYAMLConfiguration &yaml,
                    const std::vector<std::string> &property_path,
                    const std::map<std::string, T> &property_map, T &property) {
  std::string property_string;
  yaml.GetProperty(property_path, property_string, true);

  try {
    property = property_map.at(property_string);
  } catch (std::exception &exp) {
    throw std::invalid_argument(BuildErrorMessagePath(property_path, exp));
  }
}

/* Given a value of type T, returns 1 if val>0, -1 if val<0 and 0 if val ==0 */
template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
}  // namespace yaml
}  // namespace dhfe

#endif