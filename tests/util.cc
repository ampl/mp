/*
 Test utilities

 Copyright (C) 2013 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "tests/util.h"

#include "solvers/util/error.h"
#include "solvers/util/os.h"

extern "C" {
#include "solvers/asl.h"
}
#undef filename

#include <fstream>

#ifdef _WIN32
# include <direct.h>
# define chdir _chdir
#else
# include <unistd.h>
# include <sys/wait.h>
#endif

#ifndef WIFEXITED
# define WIFEXITED(status) true
# define WEXITSTATUS(status) status
#endif

std::string ReadFile(fmt::StringRef name) {
  std::string data;
  std::ifstream ifs(name.c_str());
  enum { BUFFER_SIZE = 4096 };
  char buffer[BUFFER_SIZE];
  do {
    ifs.read(buffer, BUFFER_SIZE);
    data.append(buffer, static_cast<std::string::size_type>(ifs.gcount()));
  } while (ifs);
  return data;
}

void WriteFile(fmt::StringRef name, fmt::StringRef data) {
  std::ofstream ofs(name.c_str());
  ofs.write(data.c_str(), data.size());
}

std::string FixPath(fmt::StringRef path, char sep) {
  if (sep == '/')
    return path;
  std::string fixed = path;
  std::replace(fixed.begin(), fixed.end(), '/', sep);
  return fixed;
}

StderrRedirect::StderrRedirect(const char *filename) : saved_stderr(Stderr) {
  Stderr = std::fopen(filename, "w");
}

StderrRedirect::~StderrRedirect() {
  std::fclose(Stderr);
  Stderr = saved_stderr;
}

void ChangeDirectory(fmt::StringRef path) {
  if (chdir(path.c_str()) != 0)
    ampl::ThrowError("chdir failed, error code = {}") << errno;
}

int ExecuteShellCommand(
    fmt::StringRef command, bool throw_on_nonzero_exit_code) {
#ifdef _WIN32
  std::wstring command_str((ampl::UTF8ToUTF16(command)));
  std::replace(command_str.begin(), command_str.end(), L'/', L'\\');
  int result = _wsystem(command_str.c_str());
#else
  int result = std::system(command.c_str());
#endif
  // Check if system function failed.
  if (result == -1)
    ampl::ThrowError("system failed, error code = {}") << errno;
  // Check if process hasn't exited normally.
  if (!WIFEXITED(result)) {
    ampl::ThrowError("process hasn't exited normally, error code = {}")
      << result;
  }
  // Process exited normally - check exit code.
  int exit_code = WEXITSTATUS(result);
  if (exit_code != 0 && throw_on_nonzero_exit_code)
    ampl::ThrowError("process exited with code {}") << exit_code;
  return exit_code;
}

std::string GetExecutableDir() {
  return ampl::GetExecutablePath().remove_filename().string();
}

std::string FixBinaryPath(fmt::StringRef path) {
  ampl::path exe_path = ampl::GetExecutablePath();
  exe_path.remove_filename();
  std::string config_dir = exe_path.filename().string();
  if (config_dir != "Debug" && config_dir != "Release")
    return FixPath(path.c_str());
  ampl::path p(path.c_str(), path.c_str() + path.size());
  std::string dirname = p.remove_filename().string();
  if (!dirname.empty())
    dirname += '/';
  std::string filename = p.filename().string();
  return FixPath(dirname + config_dir + "/" + filename);
}

std::vector<std::string> Split(const std::string &s, char sep) {
  std::vector<std::string> items;
  std::string::size_type start = 0, end = 0;
  while ((end = s.find(sep, start)) != std::string::npos) {
    items.push_back(s.substr(start, end - start));
    start = end + 1;
  }
  items.push_back(s.substr(start));
  return items;
}

std::string ReplaceLine(std::string s, int line_index, const char *new_line) {
  std::string::size_type start = 0;
  while (line_index-- > 0) {
    start = s.find('\n', start);
    if (start == std::string::npos)
      throw ampl::Error("invalid line index");
    ++start;
  }
  std::string::size_type end = s.find('\n', start);
  if (end == std::string::npos)
    end = s.size();
  s.replace(start, end - start, new_line);
  return s;
}
