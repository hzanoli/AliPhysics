#include "DHFeUtils.h"

std::string dhfe::yaml::BuildErrorMessagePath(
    const std::vector<std::string> &property_path, const std::exception &exp) {
  std::string combined_path;
  for (auto &x : property_path) combined_path += x + std::string(":");

  auto path = combined_path.substr(0, combined_path.length() - 1);

  std::string error_message = "Problem to read :" + path + ". \n";
  error_message += "The following error was raised: ";
  error_message += exp.what();

  return error_message;
}