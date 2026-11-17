// Stand-in solver for nl-solver-test: checks that it is called as
// "<solver> <stub> -AMPL [options]" with <stub>.nl present, and writes the
// stub it received to <stub>.ok, so that the test can compare it
// with the one NLSolver used. It does not read the model; it writes an
// empty solution file <stub>.sol, and the options it received, one per
// line, to <stub>.args. The first option selects a behavior:
//   nosol        write no solution file (exit code 0)
//   nosol_exit7  write no solution file, exit with code 7
//   exit7        write the solution file, then exit with code 7

#include <cstdio>
#include <cstring>
#include <string>

int main(int argc, char* argv[]) {
  if (argc < 3 || std::strcmp(argv[2], "-AMPL"))
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
  if (std::fclose(ok))
    return 5;
  if (auto args = std::fopen((stub + ".args").c_str(), "wb")) {
    for (int i = 3; i < argc; ++i)
      std::fprintf(args, "%s\n", argv[i]);
    std::fclose(args);
  }
  std::string mode = argc > 3 ? argv[3] : "";
  if ("nosol" == mode)
    return 0;
  if ("nosol_exit7" == mode)
    return 7;
  if (auto sol = std::fopen((stub + ".sol").c_str(), "wb"))
    std::fclose(sol);
  else
    return 6;
  return "exit7" == mode ? 7 : 0;
}
