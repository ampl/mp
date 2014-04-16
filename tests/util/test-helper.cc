// Helper program for testing GetExecutablePath and ExecuteShellCommand.

#include "solvers/util/os.h"

int main() {
  fmt::Print("{}") << ampl::GetExecutablePath().string();
  return 42;
}
