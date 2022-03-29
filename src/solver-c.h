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

#ifndef MP_SOLVER_C_H_
#define MP_SOLVER_C_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _WIN32
# define MP_API
#else
# ifdef MP_EXPORT
#  define MP_API __declspec(dllexport)
# else
#  define MP_API __declspec(dllimport)
# endif
#endif  // _WIN32

/// DEPRECATED
/// Use ampls_c_api.h and related files

/**
 * An error.
 */
typedef struct MP_Error MP_Error;

/**
 * An AMPL solver.
 */
typedef struct MP_Solver MP_Solver;

/**
 * Creates a solver object. In case of an error it returns a null
 * pointer and, if e is non-null, stores a pointer to an error object
 * in a location pointed to by e. The returned solver object should be
 * destroyed with MP_DestroySolver once it is no longer needed.
 *
 * options: A string containing solver initialization options.
 * e: A pointer to a location where to store an error, can be null
 *    in which case the error object is discarded. If it is non-null
 *    and there was an error, the error object should be destroyed
 *    with MP_DestroyError once it is no longer needed.
 */
MP_API MP_Solver *MP_CreateSolver(const char *options, MP_Error **e);

/**
 * Destroys the solver object and deallocates memory where it was stored.
 *
 * s: The solver object to destroy.
 */
MP_API void MP_DestroySolver(MP_Solver *s);

/**
 * Returns a pointer to an object that provides information about the
 * last error that occurred in the context of solver s or a null pointer
 * if there was no error. The error object is stored in the solver
 * and has the same lifetime. There is no need to call MP_DestroyError
 * on the returned object.
 *
 * s: The solver in the context of which the error has occurred.
 */
MP_API MP_Error *MP_GetLastError(MP_Solver *s);

/**
 * Destroys the error object and deallocates memory where it was stored.
 *
 * s: The error object to destroy.
 */
MP_API void MP_DestroyError(MP_Error *e);

/**
 * Returns a pointer to the error message or a null pointer if there was an
 * error. The pointer is invalidated when the error object is destroyed.
 *
 * e: The error object to query.
 */
MP_API const char *MP_GetErrorMessage(MP_Error *e);

/**
 * Returns a pointer to the option header. The option header is a text printed
 * before option descriptions.
 *
 * s: The solver object to query.
 */
MP_API const char *MP_GetOptionHeader(MP_Solver *s);

/**
 * A solver option.
 */
typedef struct MP_SolverOption MP_SolverOption;

/**
 * Solver option flags.
 */
enum {
  /**
   * If set the option contains information about possible values that can be
   * obtained via MP_GetOptionValues.
   */
  MP_OPT_HAS_VALUES = 1
};

/**
 * Solver option information.
 */
typedef struct MP_SolverOptionInfo {
  const char *name;         /**< The option name. */
  const char *description;  /**< The option description. */
  int flags;                /**< The option flags. */
  MP_SolverOption *option;
} MP_SolverOptionInfo;

/**
 * Retrieves the information about solver options. Returns the number of
 * options if succeeded, -1 otherwise.
 *
 * s: The solver object to query.
 * options: A pointer to an array of the specified size where to store
 *          the option information. Can be null in which case the function
 *          only returns the number of options.
 * size: The size of the options array. Ignored if options is null.
 */
MP_API int MP_GetSolverOptions(
    MP_Solver *s, MP_SolverOptionInfo *options, int size);

/**
 * Information about a possible option value.
 */
struct MP_OptionValueInfo {
  const char *value;        /**< The value. */
  const char *description;  /**< The value description. */
};

/**
 * Retrieves the information about possible option values.
 * Returns the number of values if succeeded, -1 otherwise.
 *
 * s: The solver object containing the option.
 * option: The option object to query.
 * values: A pointer to an array of the specified size where to store
 *         the value information. Can be null in which case the function
 *         only returns the number of values.
 * size: The size of the values array. Ignored if values is null.
 */
MP_API int MP_GetOptionValues(MP_Solver *s,
    MP_SolverOption *option, MP_OptionValueInfo *values, int size);

/**
 * Sets the value of a string option. Returns 0 if succeeded, -1 otherwise.
 *
 * s: The solver object containing the option.
 * option: The option name.
 * value: The option value.
 */
MP_API int MP_SetStrOption(MP_Solver *s, const char *option, const char *value);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // MP_SOLVER_C_H_
