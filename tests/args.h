/*
 Definition of the Args class.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef TESTS_ARGS_H_
#define TESTS_ARGS_H_

#include <cstddef>
#include <vector>

// Helper class that copies arguments to comply with the main function
// signature and avoid unwanted modification.
class Args {
 private:
  std::size_t argc_;
  std::vector<char> store_;
  std::vector<char*> argv_;
  char **pargv_;

  void Add(const char *arg);

 public:
  explicit Args(const char *arg1, const char *arg2 = 0,
      const char *arg3 = 0, const char *arg4 = 0,
      const char *arg5 = 0, const char *arg6 = 0);

  operator char **&();
};

#endif  // TESTS_ARGS_H_
