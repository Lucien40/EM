
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char const *argv[]) {
  std::string dir("tets");
  std::filesystem::create_directory(dir.c_str());
  std::ofstream fileObservables;
  fileObservables.open("/home/lucien/Documents/Em/EM/tets/blas.out");
  fileObservables << "ho";
  return 0;
}
