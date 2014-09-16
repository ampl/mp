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

#define MP_EXPORT
#include "solver-c.h"

#include "mp/solver.h"

#include <cstring>
#include <exception>
#include <memory>

extern "C" {

// Flags for MP_Error.
enum {
  DELETE_MESSAGE = 1,  // The message has to be deleted.
  DELETE_OBJECT  = 2   // The object has to be deleted.
};

struct MP_Error {
  const char *message;
  unsigned flags;
};

struct MP_Solver {
  mp::SolverPtr solver;
  MP_Error last_error;
  explicit MP_Solver(const char *options)
  : solver(mp::CreateSolver(options)) {
    last_error.message = 0;
    last_error.flags = 0;
  }
};
}  // extern "C"

namespace {

const char OUT_OF_MEMORY_MSG[] = "out of memory";
MP_Error out_of_memory = {OUT_OF_MEMORY_MSG, 0};

// Sets the error message deleting the old one if necessary.
void SetErrorMessage(MP_Error &e, const char *message) FMT_NOEXCEPT(true) {
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

void SetError(MP_Error **e, const char *message) FMT_NOEXCEPT(true) {
  if (!e) return;
  MP_Error *error = new (std::nothrow) MP_Error;
  error->flags = 0;
  if (error) {
    SetErrorMessage(*error, message);
    error->flags |= DELETE_OBJECT;
  } else {
    error = &out_of_memory;
  }
  *e = error;
}

inline void SetError(MP_Solver *s, const char *message) FMT_NOEXCEPT(true) {
  SetErrorMessage(s->last_error, message);
}
}

extern "C" {

MP_Solver *MP_CreateSolver(const char *options, MP_Error **error) {
  try {
    return new MP_Solver(options);
  } catch (const std::exception &e) {
    SetError(error, e.what());
  } catch (...) {
    SetError(error, "unknown error");
  }
  return 0;
}

void MP_DestroySolver(MP_Solver *s) {
  delete s;  // Doesn't throw.
}

MP_Error *MP_GetLastError(MP_Solver *s) {
  return s->last_error.message ? &s->last_error : 0;  // Doesn't throw.
}

void MP_DestroyError(MP_Error *e) {
  // Doesn't throw.
  if (!e) return;
  if ((e->flags & DELETE_MESSAGE) != 0)
    delete [] e->message;
  if ((e->flags & DELETE_OBJECT) != 0)
    delete e;
}

const char *MP_GetErrorMessage(MP_Error *e) {
  return e->message;  // Doesn't throw.
}

const char *MP_GetOptionHeader(MP_Solver *s) {
  try {
    return s->solver->option_header();
  } catch (const std::exception &e) {
    SetError(s, e.what());
  } catch (...) {
    SetError(s, "unknown error");
  }
  return 0;
}

int MP_GetSolverOptions(
    MP_Solver *s, MP_SolverOptionInfo *options, int size) {
  try {
    const mp::Solver &solver = *s->solver;
    int num_options = solver.num_options();
    if (!options)
      return num_options;
    int index = 0;
    for (mp::Solver::option_iterator i = solver.option_begin(),
        e = solver.option_end(); i != e && index < size; ++i, ++index) {
      const mp::SolverOption &opt = *i;
      options[index].name = opt.name();
      options[index].description = opt.description();
      options[index].flags = opt.values().size() != 0 ? MP_OPT_HAS_VALUES : 0;
      options[index].option = reinterpret_cast<MP_SolverOption*>(
          const_cast<mp::SolverOption*>(&opt));
    }
    return num_options;
  } catch (const std::exception &e) {
    SetError(s, e.what());
  } catch (...) {
    SetError(s, "unknown error");
  }
  return -1;
}

int MP_GetOptionValues(MP_Solver *s,
    MP_SolverOption *option, MP_OptionValueInfo *values, int size) {
  try {
    mp::ValueArrayRef val =
        reinterpret_cast<mp::SolverOption*>(option)->values();
    int num_values = val.size();
    if (!values)
      return num_values;
    int index = 0;
    for (mp::ValueArrayRef::iterator
        i = val.begin(), e = val.end(); i != e && index < size; ++i, ++index) {
      values[index].value = i->value;
      values[index].description = i->description;
    }
    return num_values;
  } catch (const std::exception &e) {
    SetError(s, e.what());
  } catch (...) {
    SetError(s, "unknown error");
  }
  return -1;
}

int MP_RunSolver(MP_Solver *s, int, char **argv) {
  try {
    return s->solver->Run(argv);
  } catch (const std::exception &e) {
    SetError(s, e.what());
  } catch (...) {
    SetError(s, "unknown error");
  }
  return -1;
}
}  // extern "C"
