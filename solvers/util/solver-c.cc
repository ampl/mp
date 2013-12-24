/*
 A C interface to an AMPL solver.

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

#include "solvers/util/solver-c.h"

#include "solvers/util/solver.h"

#include <cstring>
#include <exception>
#include <memory>

extern "C" {

// Flags for ASL_Error.
enum {
  DELETE_MESSAGE = 1,  // The message has to be deleted.
  DELETE_OBJECT  = 2   // The object has to be deleted.
};

struct ASL_Error {
  const char *message;
  unsigned flags;
};

struct ASL_Solver {
  ampl::SolverPtr solver;
  ASL_Error last_error;
  ASL_Solver() : solver(ampl::CreateSolver()) {}
};
}  // extern "C"

namespace {

const char OUT_OF_MEMORY_MSG[] = "out of memory";
ASL_Error out_of_memory = {OUT_OF_MEMORY_MSG};

// Sets the error message deleting the old one if necessary.
void SetErrorMessage(ASL_Error &e, const char *message) FMT_NOEXCEPT(true) {
  if (e.message && (e.flags & DELETE_MESSAGE) != 0)
    delete [] e.message;
  char *error_message = new (std::nothrow) char[std::strlen(message) + 1];
  if (error_message != 0) {
    e.flags |= DELETE_MESSAGE;
    e.message = error_message;
    std::strcpy(error_message, message);
  } else {
    e.flags &= ~DELETE_MESSAGE;
    e.message = OUT_OF_MEMORY_MSG;
  }
}

void SetError(ASL_Error **e, const char *message) FMT_NOEXCEPT(true) {
  if (!e) return;
  ASL_Error *error = new (std::nothrow) ASL_Error();
  if (error) {
    SetErrorMessage(*error, message);
    error->flags |= DELETE_OBJECT;
  } else {
    error = &out_of_memory;
  }
  *e = error;
}

inline void SetError(ASL_Solver *s, const char *message) FMT_NOEXCEPT(true) {
  SetErrorMessage(s->last_error, message);
}
}

extern "C" {

ASL_Solver *ASL_CreateSolver(ASL_Error **error) {
  try {
    return new ASL_Solver();
  } catch (const std::exception &e) {
    SetError(error, e.what());
  } catch (...) {
    SetError(error, "unknown error");
  }
  return 0;
}

void ASL_DestroySolver(ASL_Solver *s) {
  delete s;  // Doesn't throw.
}

ASL_Error *ASL_GetLastError(ASL_Solver *s) {
  return &s->last_error;  // Doesn't throw.
}

void ASL_DestroyError(ASL_Error *e) {
  // Doesn't throw.
  if (!e) return;
  if ((e->flags & DELETE_MESSAGE) != 0)
    delete [] e->message;
  if ((e->flags & DELETE_OBJECT) != 0)
    delete e;
}

const char *ASL_GetErrorMessage(ASL_Error *e) {
  return e->message;  // Doesn't throw.
}

int ASL_GetSolverOptions(ASL_Solver *s, ASL_SolverOption *options, int size) {
  try {
    const ampl::Solver &solver = *s->solver;
    int num_options = solver.num_options();
    if (!options)
      return num_options;
    int index = 0;
    for (ampl::Solver::option_iterator i = solver.option_begin(),
        e = solver.option_end(); i != e && index < size; ++i, ++index) {
      options[index].name = i->second->name();
      options[index].description = i->second->description();
    }
    return num_options;
  } catch (const std::exception &e) {
    SetError(s, e.what());
  } catch (...) {
    SetError(s, "unknown error");
  }
  return -1;
}

int ASL_RunSolver(ASL_Solver *s, int, char **argv) {
  try {
    return s->solver->Run(argv);
  } catch (const std::exception &e) {
    SetError(s, e.what());
  } catch (...) {
    SetError(s, "unknown error");
  }
  return -1;
}

// TODO: test

}  // extern "C"
