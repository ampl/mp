#include "tests/util.h"
#include "solvers/asl.h"
#undef filename

#include <cstring>
#include <fstream>

std::string ReadFile(const char *name) {
  std::string data;
  std::ifstream ifs(name);
  enum { BUFFER_SIZE = 4096 };
  char buffer[BUFFER_SIZE];
  do {
    ifs.read(buffer, BUFFER_SIZE);
    data.append(buffer, static_cast<std::string::size_type>(ifs.gcount()));
  } while (ifs);
  return data;
}

void WriteFile(const char *name, const char *data) {
  std::ofstream ofs(name);
  ofs.write(data, std::strlen(data));
}

StderrRedirect::StderrRedirect(const char *filename) : saved_stderr(Stderr) {
  Stderr = std::fopen(filename, "w");
}

StderrRedirect::~StderrRedirect() {
  std::fclose(Stderr);
  Stderr = saved_stderr;
}
