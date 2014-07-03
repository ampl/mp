// Helper program for testing GetExecutablePath and ExecuteShellCommand.

#include "solvers/util/os.h"

int main() {
  fmt::print("{}", ampl::GetExecutablePath().string());
  return 42;
}
