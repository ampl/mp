/*
 SOL reader2.

 SOL is a format for representing solutions of optimization problems
 such as linear, quadratic, nonlinear, complementarity
 and constraint programming problems in discrete or continuous variables,
 used by AMPL. It is not documented yet.

 Usage:
   MySOLHandler handler;
	 mp::NLUtils utils;

	 auto status = ReadSOLFile(name, handler, utils);
	 if (status.first) {
     printf("Error: %s\n", status.second.c_str());
		 exit(EXIT_FAILURE);
   }

 where handler is an object that receives information on solution
 components.

 Copyright (C) 2023 AMPL Optimization Inc.

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

 Author: Gleb Belov
 */

#ifndef SOLREADER2_H
#define SOLREADER2_H

#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <vector>

#include "mp/nl-solver-basics-c.h"
#include "mp/sol-handler.h"
#include "mp/nl-utils.h"

namespace mp {

/// Vector reader
template <class Value>
class VecReader {
public:
  /// Construct
  VecReader(FILE* f, int b, int n)
    : f_(f), binary_(b), n_(n) { }

  /// Destruct
  ~VecReader() { assert(!n_); }

  /// Typedef value_type
  using value_type = Value;

  /// Size remaining to read
  int Size() const { return n_; }

  /// Read next value
  Value ReadNext();

  /// Set error status
  void SetError(NLW2_SOLReadResultCode res, std::string msg)
  { rr_ = res; err_msg_ = std::move(msg); n_ = 0; }

  /// status
  NLW2_SOLReadResultCode ReadResult() const { return rr_; }

  /// error message
  const std::string& ErrorMessage() const
  { return err_msg_; }


private:
	FILE* f_;
	int binary_;
  int n_;

  NLW2_SOLReadResultCode rr_{NLW2_SOLRead_OK};
  std::string err_msg_;
};


/// Suffix info
class SuffixInfo {
public:
  /// Construct
  SuffixInfo(int k, std::string n,
             std::string t)
    : kind_(k), name_(std::move(n)),
      table_(std::move(t)) { }
  /// Suffix kind
  int Kind() const { return kind_; }
  /// Suffix name
  const std::string& Name() const { return name_; }
  /// Suffix table
  const std::string& Table() const { return table_; }

private:
  int kind_;
  std::string name_, table_;
};


/// Typedef SparseVecReader
template <class Element>
using SparseVecReader
    = VecReader< std::pair<int, Element> >;


/// Class SuffixReader<>
template <class Element>
class SuffixReader
    : public SparseVecReader<Element> {
public:
  /// Construct
  SuffixReader(SuffixInfo si, FILE* f, int b, int n)
    : SparseVecReader<Element>(f, b, n),
      si_(std::move(si)) { }
  /// Suffix info
  const SuffixInfo& SufInfo() const { return si_; }

private:
  SuffixInfo si_;
};

/// Declare struct SufRead
struct SufRead;


/// SOL reader.
///
/// AMPL SOL is a format for representing solutions of optimization problems
/// such as linear, quadratic, nonlinear, complementarity
/// and constraint programming problems in discrete or continuous variables,
/// used by AMPL. It is not documented yet.
///
/// Usage: recommended via mp::NLSOL class.
/// If you need to read solutions separately, proceed as follows:
///
///   MySOLHandler handler;
///   mp::NLUtils utils;
///
///   auto status = ReadSOLFile(name, handler, utils);
///   if (status.first) {
///     printf("Error: %s\n", status.second.c_str());
///     exit(EXIT_FAILURE);
///   }
///
/// @param  SOLHandler: class that receives information on solution
/// components.
template <typename SOLHandler>
class SOLReader2 {
public:
  /// Constructor.
	SOLReader2(SOLHandler& sh, NLUtils& utl) :
			solh_(sh), utils_(utl), hdr_(sh.Header()) { }

  /// Read .sol file.
  NLW2_SOLReadResultCode ReadSOLFile(const std::string& name);

  /// Error message, if any.
  /// Warnings are printed independently
  /// via NLUtils::log_warning().
  const std::string& ErrorMessage(NLW2_SOLReadResultCode ) const
  { return err_msg_; }

  /// Internal AMPL result code
  int internal_rv() const { return internal_rv_; }


protected:
  /// Get const handler&
  const SOLHandler& Handler() const { return solh_; }
  /// Get handler&
  SOLHandler& Handler() { return solh_; }

  /// Get const utils&
  const NLUtils& Utils() const { return utils_; }
  /// Get utils&
  NLUtils& Utils() { return utils_; }

	/// Get const NLHeader&
	const NLHeader& Header() const { return hdr_; }

	/// NumVars
	int NumVars() const { return Header().num_vars; }

	/// NumAlgCons
	int NumAlgCons() const { return Header().num_algebraic_cons; }

	/// NumLogCons
	int NumLogCons() const { return Header().num_logical_cons; }

  /// Check result of a [vector] reader
  template <class Reader>
  bool CheckReader(const Reader& rd, NLW2_SOLReadResultCode& rr) {
    if (NLW2_SOLRead_Early_EOF == rd.ReadResult())
      return (rr=ReportEarlyEof(), false);
    if (NLW2_SOLRead_Bad_Line == rd.ReadResult())
      return (rr=ReportBadLine(rd.ErrorMessage()), false);
    if (rd.Size()) {       // unfinished
      serror("vector not read completely");
      return (rr=ReportBadFormat(), false);
    }
    if (NLW2_SOLRead_OK != rd.ReadResult()) {
      serror( rd.ErrorMessage().c_str() );
      ReportBadFormat();
      return (rr=rd.ReadResult(), false);
    }
    return true;
  }

  /// Read suffixes, binary
  NLW2_SOLReadResultCode bsufread(FILE* f);

  /// Read suffixes, text
  NLW2_SOLReadResultCode gsufread(FILE* f);

  /// Suffix head check
  int sufheadcheck(SufRead *sr);

  /// Report Early Eof
  NLW2_SOLReadResultCode ReportEarlyEof();

  /// Report bad format
  NLW2_SOLReadResultCode ReportBadFormat();

  /// Report bad line
  NLW2_SOLReadResultCode ReportBadLine(const std::string& line);

  /// Save error message
  void serror(const char* format, ...);


#ifndef Long
  using Long = int;
  using uLong = unsigned int;
  /** cfront screws up when we use typedef Long uiolen; */
  /** It renders Long as long rather than Long. */
  using uiolen = uLong;
#else
  static_assert (sizeof(Long)==sizeof(int), "Unusable Long");
#endif

  using real = double;

#ifdef CRAY
  typedef int integer;
#else
  typedef uLong integer;
#endif

private:
  SOLHandler& solh_;
  NLUtils& utils_;
	NLHeader hdr_;

  std::string solve_msg_;
  std::string err_msg_;
  int internal_rv_ {0};
  NLW2_SOLReadResultCode readresult_ {NLW2_SOLRead_Result_Not_Set};

  const char* stub_ {nullptr};
  const char *bkind[2] = { "ASCII", "binary" };
  char *b1, buf[512], *from_solve, *s, *se;
  Long Objno[2] = {-2, -2};	/* Objno[1] = solve_result_num */
  Long nOpts, Lt, Options[14], *z;
  uiolen L, L1, L2;
  int binary, bs, c, have_options, hs, i, j, j1, je;
  int need_vbtol, newsufs, nsv, objno=-2,
  presolved, rv, slin, solmsg;
  double vbtol, x;
  size_t len, n, n1, nbs, ui;
  int ndebug, sstatus_seen;
};


/// Read SOL file.
template <class SOLHandler>
inline std::pair<NLW2_SOLReadResultCode, std::string>
ReadSOLFile(
    const std::string& name,
		SOLHandler& solh, NLUtils& utl, int* p_internal_rv=nullptr) {
  SOLReader2<SOLHandler> solr(solh, utl);
	auto res = solr.ReadSOLFile(name);
	if (p_internal_rv)
		*p_internal_rv = solr.internal_rv();
	return { res, solr.ErrorMessage(res) };
}


}  // namespace mp

#endif // SOLREADER2_H
