/**
 Common definitions for NL reader and NL writer.

 Copyright (C) 2014 - 2023 AMPL Optimization Inc.

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

 Authors: Victor Zverovich, Gleb Belov
 */

#ifndef NLHEADER_H
#define NLHEADER_H

#include <cstddef>
#include <algorithm>
#include <array>
#include <string>

/// Namespace mp.
namespace mp {


/** Information about an optimization problem.
 *  Full documentation on the NL format:
 *  technical report "Writing .nl Files"
 *  (https://ampl.github.io/nlwrite.pdf.)
 */
struct ProblemInfo {
	/** Total number of variables. */
	int num_vars {0};

	/**
		Number of algebraic constraints including ranges and equality constraints.
		It doesn't include logical constraints.
	 */
	int num_algebraic_cons {0};

	/** Total number of objectives. */
	int num_objs {0};

	/** Number of ranges (constraints with -Infinity < LHS < RHS < Infinity). */
	int num_ranges {0};

	/**
		Number of equality constraints or -1 if unknown (AMPL prior to 19970627).
	 */
	int num_eqns {0};

	/** Number of logical constraints. */
	int num_logical_cons {0};

	/** Number of random variables. */
	int num_rand_vars {0};

	/** Number of random defined variables. */
	int num_rand_common_exprs {0};

	/** Number of random constraints. */
	int num_rand_cons {0};

	/** Number of random objectives. */
	int num_rand_objs {0};

	/** Number of random calls. */
	int num_rand_calls {0};

	/** Number of stages. */
	int num_stages {0};

  /** Returns the number of integer variables (includes binary variables). */
	int num_integer_vars() const {
		return num_linear_binary_vars + num_linear_integer_vars +
				num_nl_integer_vars_in_both + num_nl_integer_vars_in_cons +
				num_nl_integer_vars_in_objs;
	}

	/** Returns the number of continuous variables. */
	int num_continuous_vars() const { return num_vars - num_integer_vars(); }

	// Nonlinear and complementarity information
	// -----------------------------------------

	/** Total number of nonlinear constraints. */
	int num_nl_cons {0};

	/** Total number of nonlinear objectives. */
	int num_nl_objs {0};

	/** Total number of complementarity conditions. */
	int num_compl_conds {0};

	/** Number of nonlinear complementarity conditions. */
	int num_nl_compl_conds {0};

	/** Number of complementarities involving double inequalities. */
	int num_compl_dbl_ineqs {0};

	/** Number of complemented variables with a nonzero lower bound. */
	int num_compl_vars_with_nz_lb {0};

	// Information about network constraints
	// -------------------------------------

	/** Number of nonlinear network constraints. */
	int num_nl_net_cons {0};

	/** Number of linear network constraints. */
	int num_linear_net_cons {0};

	// Information about nonlinear variables
	// -------------------------------------

	/**
		Number of nonlinear variables in constraints including nonlinear
		variables in both constraints and objectives.
	 */
	int num_nl_vars_in_cons {0};

	/**
		Number of nonlinear variables in objectives including nonlinear
		variables in both constraints and objectives.
	 */
	int num_nl_vars_in_objs {0};

	/** Number of nonlinear variables in both constraints and objectives. */
	int num_nl_vars_in_both {0};

	// Miscellaneous
	// -------------

	/** Number of linear network variables (arcs). */
	int num_linear_net_vars {0};

	/** Number of functions. */
	int num_funcs {0};

	// Information about discrete variables
	// ------------------------------------

	/** Number of linear binary variables. */
	int num_linear_binary_vars {0};

	/** Number of linear non-binary integer variables. */
	int num_linear_integer_vars {0};

	/**
		Number of integer nonlinear variables in both constraints and objectives.
	 */
	int num_nl_integer_vars_in_both {0};

	/** Number of integer nonlinear variables just in constraints. */
	int num_nl_integer_vars_in_cons {0};

	/** Number of integer nonlinear variables just in objectives. */
	int num_nl_integer_vars_in_objs {0};

	// Information about nonzeros
	// --------------------------

	/** Number of nonzeros in constraints' Jacobian. */
	std::size_t num_con_nonzeros {0};

	/** Number of nonzeros in all objective gradients. */
	std::size_t num_obj_nonzeros {0};

  // Information about names.
  // Does not have to be filled for NLWriter2.
	// -----------------------

  /** Length of longest constraint or objective name
   *  if names are available. */
	int max_con_name_len {0};

	/** Length of longest variable name if names are available. */
	int max_var_name_len {0};

	// Information about common expressions
	// ------------------------------------

	/**
		Number of common expressions that appear both in constraints
		and objectives.
	 */
	int num_common_exprs_in_both {0};

	/**
		Number of common expressions that appear in multiple constraints
		and don't appear in objectives.
	 */
	int num_common_exprs_in_cons {0};

	/**
    Number of common expressions that appear in multiple objectives
		and don't appear in constraints.
	 */
	int num_common_exprs_in_objs {0};

	/**
		Number of common expressions that only appear in a single constraint
		and don't appear in objectives.
	 */
	int num_common_exprs_in_single_cons {0};

	/**
		Number of common expressions that only appear in a single objective
		and don't appear in constraints.
	 */
	int num_common_exprs_in_single_objs {0};

	/** Returns the total number of common expressions. */
	int num_common_exprs() const {
		return num_common_exprs_in_both + num_common_exprs_in_cons +
				num_common_exprs_in_objs + num_common_exprs_in_single_cons +
				num_common_exprs_in_single_objs;
	}
};


enum {
	/** Maximum number of options reserved for AMPL use in NL and SOL formats. */
	MAX_AMPL_OPTIONS = 9
};

enum {
	VBTOL_OPTION_INDEX = 1,
	USE_VBTOL_FLAG     = 3
};

/// Namespace arith
namespace arith {

/** Floating-point arithmetic kind. */
enum Kind {

	/** Unknown floating-point arithmetic. */
	UNKNOWN = 0,

	/**
		\rst
		Standard `IEEE-754 floating point
		<http://en.wikipedia.org/wiki/IEEE_floating_point>`_ - little endian.
		\endrst
	 */
	IEEE_LITTLE_ENDIAN = 1,

	/** Standard IEEE-754 floating point - big endian. */
	IEEE_BIG_ENDIAN = 2,

	/**
		\rst
		`IBM floating point
		<http://en.wikipedia.org/wiki/IBM_Floating_Point_Architecture>`_.
		\endrst
	 */
	IBM = 3,

	/** VAX floating point (legacy). */
	VAX = 4,

	/** Cray floating point. */
	CRAY = 5,

	/** Last floating point. */
	LAST = CRAY
};

/// Returns floating-point arithmetic kind used on the current system.
Kind GetKind();

/// IsIEEE
inline bool IsIEEE(arith::Kind k) {
	return k == IEEE_LITTLE_ENDIAN || k == IEEE_BIG_ENDIAN;
}
}  // namespace arith


/**
	\rst
	An NL `header <http://en.wikipedia.org/wiki/Header_(computing)>`_
	which contains information about problem dimensions, such as the number of
	variables and constraints, and the input format.

	Base class: `mp::ProblemInfo`
	\endrst
 */
struct NLHeader : ProblemInfo {
	/** Input/output format */
	enum Format {
		/**
			Text format. The text format is fully portable meaning that an .nl file
			can be written on a machine of one architecture and then read on a
			machine of a different architecture.
		 */
		TEXT = 0,

		/**
			Binary format. The binary format is not generally portable and should
			normally be used on a single machine.
		 */
		BINARY = 1
	};

	/** Input/output format. */
	Format format {TEXT};

	/** The number of options reserved for AMPL use. */
	int num_ampl_options {3};

	/**
    Values of options reserved for AMPL use.
    Leave the default values if not using AMPL.
	 */
  long ampl_options[MAX_AMPL_OPTIONS] = {1, 1, 0};

	/**
    Extra info for writing a solution reserved for AMPL use.
    Leave the default value if not using AMPL.
	 */
	double ampl_vbtol {0.0};

	/**
	 * Problem name.
	 */
	std::string prob_name;

	/**
		\rst
		Floating-point arithmetic kind used with binary format to check
		if an .nl file is written using a compatible representation of
		floating-point numbers. It is not used with the text format and normally
		set to `mp::arith::UNKNOWN` there.
		\endrst
	 */
  arith::Kind arith_kind {arith::IEEE_LITTLE_ENDIAN};

	/** Flags. */
	enum {
		/** Flag that specifies whether to write output suffixes to a .sol file. */
		WANT_OUTPUT_SUFFIXES = 1
	};

	/**
		\rst
		Flags. Can be either 0 or `mp::NLHeader::WANT_OUTPUT_SUFFIXES`.
		\endrst
	 */
  int flags {WANT_OUTPUT_SUFFIXES};

	NLHeader()
    : ProblemInfo(),
      format(TEXT),
      num_ampl_options(3), ampl_vbtol(0),
      arith_kind(arith::IEEE_LITTLE_ENDIAN),
      flags(WANT_OUTPUT_SUFFIXES) {
		std::fill(ampl_options, ampl_options + MAX_AMPL_OPTIONS, 0);
		std::array<long, 3> opt_default {1, 1, 0};
		std::copy(opt_default.begin(), opt_default.end(), ampl_options);
	}
};

}  // namespace mp

#endif // NLHEADER_H
