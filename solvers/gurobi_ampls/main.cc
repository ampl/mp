/// This simple driver uses GurobiAMPLS C interface

#include <cstdio>

extern "C" {
  #include "gurobi-ampls.h"
}

int main(int argc, char *argv[])
{
  if (argc<2) {
    std::printf("Usage: gurobi-ampls model.nl [\"option1=val ...\"]  \n"
                " \n"
                "Include 'wantsol=1' for a .sol file. \n");
    return -1;
  }

  const char* slv_opt = (argc<3) ? nullptr : argv[2];
  return RunGurobiAMPLS(argv[1], slv_opt);
}

