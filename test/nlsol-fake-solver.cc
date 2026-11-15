// Stand-in solver for nl-solver-test: checks that it is called as
// "<solver> <stub> -AMPL" with <stub>.nl present, and writes the
// stub it received to <stub>.ok, so that the test can compare it
// with the one NLSolver used. It does not read the model or
// write a solution.

#include <cstdio>
#include <cstring>
#include <string>

int main(int argc, char* argv[]) {
  if (argc != 3 || std::strcmp(argv[2], "-AMPL"))
    return 2;
  std::string stub = argv[1];
  if (auto nl = std::fopen((stub + ".nl").c_str(), "rb"))
    std::fclose(nl);
  else
    return 3;
  auto ok = std::fopen((stub + ".ok").c_str(), "wb");
  if (!ok)
    return 4;
  std::fputs(stub.c_str(), ok);
  return std::fclose(ok) ? 5 : 0;
}
