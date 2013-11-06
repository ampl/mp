// This program prints executable path for testing purposes.

#include "solvers/util/os.h"

int main() {
  fmt::Print("{}") << ampl::GetExecutablePath().string();
}
