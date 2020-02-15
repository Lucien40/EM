
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char const *argv[]) {
  // raw-string literal example with the literal
  // made up of separate, concatenated literals
  std::string s =
      R"(abc)"
      R"( followed by not a newline: \n)"
      " which is then followed by a non-raw literal that's concatenated \n with"
      " an embedded non-raw newline";

  std::cout << s << std::endl;

  return 0;
}
