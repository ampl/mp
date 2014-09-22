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

#include "util.h"

#include "mp/error.h"
#include "mp/os.h"

#include <cerrno>
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
  fmt::BufferedFile file(name, "wb");
  std::size_t size = data.size();
  if (std::fwrite(data.c_str(), 1, size, file.get()) != size)
    throw fmt::SystemError(errno, "cannot write file {}", name);
}

std::string FixPath(fmt::StringRef path, char sep) {
  if (sep == '/')
    return path;
  std::string fixed = path;
  std::replace(fixed.begin(), fixed.end(), '/', sep);
  return fixed;
}

void ChangeDirectory(fmt::StringRef path) {
  if (chdir(path.c_str()) != 0)
    throw mp::Error("chdir failed, error code = {}", errno);
}

int ExecuteShellCommand(
    fmt::StringRef command, bool throw_on_nonzero_exit_code) {
#ifdef _WIN32
  std::wstring command_str = fmt::internal::UTF8ToUTF16(command).str();
  std::replace(command_str.begin(), command_str.end(), L'/', L'\\');
  int result = _wsystem(command_str.c_str());
#else
  int result = std::system(command.c_str());
#endif
  // Check if system function failed.
  if (result == -1)
    throw mp::Error("system failed, error code = {}", errno);
  // Check if process hasn't exited normally.
  if (!WIFEXITED(result)) {
    throw mp::Error(
        "process hasn't exited normally, error code = {}", result);
  }
  // Process exited normally - check exit code.
  int exit_code = WEXITSTATUS(result);
  if (exit_code != 0 && throw_on_nonzero_exit_code)
    throw mp::Error("process exited with code {}", exit_code);
  return exit_code;
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
      throw mp::Error("invalid line index");
    ++start;
  }
  std::string::size_type end = s.find('\n', start);
  if (end == std::string::npos)
    end = s.size();
  s.replace(start, end - start, new_line);
  return s;
}

mp::NLHeader MakeTestHeader() {
  mp::NLHeader h = mp::NLHeader();
  h.format = mp::NLHeader::TEXT;
  h.num_options = 9;
  int options[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
  for (int i = 0; i < h.num_options; ++i)
    h.options[i] = options[i];
  h.ampl_vbtol = 1.23;

  h.num_vars = 29;
  h.num_algebraic_cons = 47;
  h.num_objs = 37;
  h.num_ranges = 41;
  h.num_eqns = 43;
  h.num_logical_cons = 31;

  h.num_nl_cons = 53;
  h.num_nl_objs = 59;
  h.num_compl_conds = 67;
  h.num_nl_compl_conds = 61;
  h.num_compl_dbl_ineqs = 71;
  h.num_compl_vars_with_nz_lb = 73;

  h.num_nl_net_cons = 79;
  h.num_linear_net_cons = 83;

  h.num_nl_vars_in_cons = 89;
  h.num_nl_vars_in_objs = 97;
  h.num_nl_vars_in_both = 101;

  h.num_linear_net_vars = 103;
  h.num_funcs = 107;
  h.arith_kind = mp::arith::IEEE_LITTLE_ENDIAN;
  h.flags = 109;

  h.num_linear_binary_vars = 113;
  h.num_linear_integer_vars = 127;
  h.num_nl_integer_vars_in_both = 131;
  h.num_nl_integer_vars_in_cons = 137;
  h.num_nl_integer_vars_in_objs = 139;

  h.num_con_nonzeros = 149;
  h.num_obj_nonzeros = 151;

  h.max_con_name_len = 157;
  h.max_var_name_len = 163;

  h.num_common_exprs_in_both = 167;
  h.num_common_exprs_in_cons = 173;
  h.num_common_exprs_in_objs = 179;
  h.num_common_exprs_in_single_cons = 181;
  h.num_common_exprs_in_single_objs = 191;
  return h;
}
