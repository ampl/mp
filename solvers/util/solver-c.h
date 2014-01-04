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

#ifndef ASL_SOLVER_C_H_
#define ASL_SOLVER_C_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
 * An error.
 */
typedef struct ASL_Error ASL_Error;

/**
 * An AMPL solver.
 */
typedef struct ASL_Solver ASL_Solver;

/**
 * Creates a solver object. In case of an error it returns a null
 * pointer and, if e is non-null, stores a pointer to an error object
 * in a location pointed to by e. The returned solver object should be
 * destroyed with ASL_DestroySolver once it is no longer needed.
 *
 * e: A pointer to a location where to store an error, can be null
 *    in which case the error object is discarded. If it is non-null
 *    and there was an error, the error object should be destroyed
 *    with ASL_DestroyError once it is no longer needed.
 */
ASL_Solver *ASL_CreateSolver(ASL_Error **e);

/**
 * Destroys the solver object and deallocates memory where it was stored.
 *
 * s: The solver object to destroy.
 */
void ASL_DestroySolver(ASL_Solver *s);

/**
 * Returns a pointer to an object that provides information about the
 * last error that occurred in the context of solver s or a null pointer
 * if there was no error. The error object is stored in the solver
 * and has the same lifetime. There is no need to call ASL_DestroyError
 * on the returned object.
 *
 * s: The solver in the context of which the error has occurred.
 */
ASL_Error *ASL_GetLastError(ASL_Solver *s);

/**
 * Destroys the error object and deallocates memory where it was stored.
 *
 * s: The error object to destroy.
 */
void ASL_DestroyError(ASL_Error *e);

/**
 * Returns a pointer to the error message. The pointer is invalidated
 * when the error object is destroyed.
 *
 * e: The error object to query.
 */
const char *ASL_GetErrorMessage(ASL_Error *e);


typedef struct ASL_SolverOption ASL_SolverOption;

/**
 * Solver option information.
 */
typedef struct ASL_SolverOptionInfo {
  const char *name;         /**< The option name. */
  const char *description;  /**< The option description. */
  ASL_SolverOption *option;
} ASL_SolverOptionInfo;

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
int ASL_GetSolverOptions(
    ASL_Solver *s, ASL_SolverOptionInfo *options, int size);

/**
 * A value of an enumerated solver option.
 */
struct ASL_EnumOptionValue {
  const char *value;        /**< The value. */
  const char *description;  /**< The value description. */
};

/**
 * Retrieves the information about possible values of an enumerated option.
 * Returns the number of values if succeeded, -1 otherwise.
 *
 * s: The solver object containing the option.
 * s: The option object to query.
 * values: A pointer to an array of the specified size where to store
 *         the value information. Can be null in which case the function
 *         only returns the number of values.
 * size: The size of the values array. Ignored if values is null.
 */
int ASL_GetOptionValues(ASL_Solver *s,
    ASL_SolverOption *option, ASL_EnumOptionValue *values, int size);

/**
 * Runs the solver. This is a programmatic alternative to running the solver
 * executable. It processes command-line arguments and solver options from
 * argv and environment variables and, if the arguments contain the filename
 * stub, reads the problem, solves it and writes the solution. Returns the
 * solver exit code.
 *
 * s: The solver to run.
 * argc: The number of arguments.
 * argv: The array of arguments of size argc + 1 with argv[argc] being
 *       a null pointer.
 */
int ASL_RunSolver(ASL_Solver *s, int argc, char **argv);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // ASL_SOLVER_C_H_
