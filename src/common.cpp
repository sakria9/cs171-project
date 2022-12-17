#include "common.h"

std::string GetFilePath(const std::string& target, int depth/* = 5*/) {
  std::string path = target;
  for (int i = 0; i < depth; ++i) {
    FILE *file = fopen(path.c_str(), "r");
    if (file) {
      fclose(file);
      return path;
    }
    path = "../" + path;
  }
  std::cerr << "GetFilePath(): Failed to get file: " + target << std::endl;
  abort();
}
