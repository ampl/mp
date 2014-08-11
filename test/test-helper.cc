// Helper program for testing GetExecutablePath and ExecuteShellCommand.

#include "mp/os.h"

int main() {
  fmt::print("{}\n", mp::GetExecutablePath().string());
  return 42;
}
