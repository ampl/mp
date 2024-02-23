/*
 Copyright (C) 2024 AMPL Optimization Inc.

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
#ifndef SOLREADER2_HPP
#define SOLREADER2_HPP

#include <cstdio>

#include "mp/sol-reader2.h"

namespace mp {

/// Parse double
inline int decstring(const char *buf, double *val) {
  char *be;
  int c;
	*val = std::strtod(buf, &be);
  return be <= buf
			|| (((c = be[-1]) < '0' || c > '9') && c != '.');
}

/// Read a double from text or binary file
inline NLW2_SOLReadResultCode Read(
    FILE* f, int binary, double& v, std::string& err) {
  err.resize(512);
  if (binary ? !std::fread((char *)&v, sizeof(double), 1, f)
       : !fgets((char*)err.data(), err.size()-1, f))
    return NLW2_SOLRead_Early_EOF;
  if (!binary && decstring(err.data(), &v))
    return NLW2_SOLRead_Bad_Line;
  return NLW2_SOLRead_OK;
}


/// Read a pair<int, El> from text or binary file
template <class El>
inline NLW2_SOLReadResultCode Read(
    FILE* f, int binary,
    std::pair<int, El>& v, std::string& err) {
  err.resize(512);
  if (binary) {
    if (fread(&v.first, sizeof(int), 1, f) != 1
        || fread(&v.second, sizeof(El), 1, f) != 1)
      return NLW2_SOLRead_Early_EOF;
  } else {
    if (!fgets((char*)err.data(), err.size()-1, f))
      return NLW2_SOLRead_Early_EOF;
    const char *s;
    char *se;
    v.first = (int)strtol(s = err.c_str(), &se, 10);
    if (se <= s)
      return NLW2_SOLRead_Bad_Line;
    auto el = strtod(s = se, &se);
    if (se <= s)
      return NLW2_SOLRead_Bad_Line;
    v.second = (El)el;
  }
  return NLW2_SOLRead_OK;
}



template <class Value>
Value VecReader<Value>::ReadNext() {
  Value v;
  --n_;                // decrement counter
  assert(n_ >= 0);
  if (NLW2_SOLRead_OK !=
      (rr_ = Read(f_, binary_, v, err_msg_))) {
    n_ = 0;            // return error status
  }
  return v;
}


/// Some refurbished old code.
template <class SOLHandler>
NLW2_SOLReadResultCode
SOLReader2<SOLHandler>::ReadSOLFile(
    const std::string& name) {
  File file;
  stub_ = name.c_str();

  file.Open(stub_, "rb");
  internal_rv_ = 997;
  if (!file) {
    serror("can't open '%s'", stub_);
    internal_rv_ = 998;
    return NLW2_SOLRead_Fail_Open;
  }
  FILE* f = file.GetHandle();
  if (fread((char *)&L, sizeof(uiolen), 1, f)
      && L == 6) {
    /* binary files may be written by Fortran unformatted writes */
    binary = 1;
    if (!fread(buf, 6, 1, f)
        || strncmp(buf,"binary",6)
        || !fread((char *)&L, sizeof(uiolen), 1, f)
        || L != 6) {
      return ReportBadFormat();
    }
  }
  else {
    binary = 0;
    rewind(f);
  }
  have_options = need_vbtol = 0;
  bs = 1;	/* omit initial backspaces */
  nbs = 0; /* number of backspaces omitted */
  if (binary) {     /////////////// BINARY FORMAT //////////
    for(;;) {       /////////////// SOLVE MESSAGE //////////
      if (!fread((char *)&L,sizeof(uiolen),1,f))
        return ReportEarlyEof();
      if ((L1 = L)) {
        do {
          n = L < sizeof(buf) ? (int)L : (int)sizeof(buf);
          L -= n;
          if (!fread(buf, n, 1, f))
            return ReportEarlyEof();
          if (!L) {
            while(n > 0) {
              if (buf[--n] != ' ') {
                ++n;
                break;
              }
            }
          }
          b1 = buf;
          n1 = n;
          if (buf[0] == '\b' && bs) {
            while(--n1 > 0
                  && *++b1 == '\b');
            nbs += n - n1;
            if (n1 > 0)
              bs = 0;
            else
              continue;
          }
          solve_msg_.append(b1, n1);
        }
        while(L);
      }
      if (!fread((char *)&L, sizeof(uiolen), 1, f))
        return ReportEarlyEof();
      if (L != L1)
        return ReportBadFormat();
      if (!L)
        break;
    }
		L1 = NumAlgCons() * sizeof(real);
    if (!fread((char *)&L, sizeof(uiolen), 1, f))
      return ReportEarlyEof();
    L2 = L - (8*sizeof(Long) + 7);
    if (L2 <= 4*sizeof(Long) + sizeof(real)
        && !(L2 & (sizeof(Long)-1))) {
      /////////////// check for Options ///////////////
      if (!fread(buf, 7, 1, f))
        return ReportEarlyEof();
      if (strncmp(buf, "Options", 7))
        return ReportBadFormat();
      if (!fread((char *)Options, sizeof(Long), 4, f))
        return ReportEarlyEof();
      nOpts = Options[0];
      if (nOpts < 3 || nOpts > 9) {
bad_nOpts:
        serror(
              "expected nOpts between 3 and 9; got %d: ",
              nOpts);
        return ReportBadFormat();
      }
      if (Options[2] == 3) {
        nOpts -= 2;
        need_vbtol = 1;
      }
      if (!fread((char *)(Options+4), sizeof(Long),
                 (size_t)(nOpts+1), f))
        return ReportEarlyEof();
      if (need_vbtol
          && !fread((char *)&vbtol, sizeof(real), 1, f))
        return ReportEarlyEof();
      if (!fread((char *)&L2, sizeof(uiolen), 1, f)
          || L != L2)
        return ReportBadFormat();
      have_options = 1;
    }
    else if (L != L1)
      return ReportBadFormat();
  }
  else {         ///////////////// TEXT FORMAT ///////////////
    for(;;) {    ///////////////// SOLVE MESSAGE /////////////
      if (!fgets(buf, sizeof(buf), f)) {
        return ReportEarlyEof();
      }
      for(se = buf; *se; se++)
        if (*se == '\r' && se[1] == '\n') {
          *se = '\n';
          *++se = 0;
          break;
        }
      if (*buf == '\n')
        break;
      n1 = se - buf;
      b1 = buf;
      if (buf[0] == '\b' && bs) {
        n = n1;
        while(--n1 > 0 && *++b1 == '\b');
        nbs += n - n1;
        if (n1 > 0)
          bs = 0;
        else
          continue;
      }
      solve_msg_.append(b1, n1);
    }
    while((j = getc(f)) == '\n' || j == '\r');
    if (j != 'O')     ////////// Check for Options //////////
      ungetc(j,f);
    else {
      if (!fgets(buf, sizeof(buf), f))
        return ReportEarlyEof();
      if (!strncmp(buf, "ptions", 6)) {
        have_options = 1;
        for(j = 0; j < 4; j++) {
          if (!fgets(buf, sizeof(buf), f))
            return ReportEarlyEof();
          Options[j] = (int)strtol(buf,&se,10);
          if (se == buf)
            return ReportBadLine(buf);
        }
        nOpts = Options[0];
        if (nOpts < 3 || nOpts > 9)
          goto bad_nOpts;
        if (Options[2] == 3) {
          nOpts -= 2;
          need_vbtol = 1;
        }
        je = (int)(nOpts+5);
        for(j = 4; j < je; j++) {
          if (!fgets(buf, sizeof(buf), f))
            return ReportEarlyEof();
          Options[j] = (int)strtol(buf,&se,10);
          if (se == buf)
            return ReportBadLine(buf);
        }
        if (need_vbtol) {
          if (!fgets(buf, sizeof(buf), f))
            return ReportEarlyEof();
          vbtol = strtod(buf,&se);
          if (se == buf)
            return ReportBadLine(buf);
        }
      }
    }
  }

  /* send termination msg */
  if (nbs) {
    // Should not have \b here but...
    auto b=solve_msg_.begin();
    while (solve_msg_.end()!=b && '\b'==*b)
      ++b;
    solve_msg_.erase(solve_msg_.begin(), b);
  }
  if (solve_msg_.size()) {
    if (binary)
      solve_msg_ += "\n";
    Handler().OnSolveMessage(
          solve_msg_.c_str(), nbs);
  }

  if (have_options) {
    z = Options + nOpts + 1;
    internal_rv_ = 996;
    // Send options to Handler:
    typename SOLHandler::AMPLOptions ao;
    ao.options_.assign(Options, Options+nOpts+5);
    ao.has_vbtol_ = need_vbtol;
    ao.vbtol_ = vbtol;
    // Handler says we should stop:
    if (auto rv = Handler().OnAMPLOptions(ao)) {
      internal_rv_ = rv;
      return NLW2_SOLRead_Bad_Options;
    }

    // Some checks.
    j = (int)z[3];
    internal_rv_ = 997;
		if (j > NumVars() || j < 0) {
      serror("Wrong NumVars %d, expected 0..%d: ",
						 j, NumVars());
      return ReportBadFormat();
      }
    j = (int)z[1];
		if (j > NumAlgCons() || j < 0) {
      serror("Wrong NumAlgCons %d, expected 0..%d: ",
						 j, NumAlgCons());
      return ReportBadFormat();
      }

    if (binary) {      // read on for binary
      L1 = j * sizeof(real);
      if (!fread((char *)&L, sizeof(uiolen), 1, f))
        return ReportEarlyEof();
      if (L != L1)
        return ReportBadFormat();
    }
  }
	else j = NumAlgCons();

  // Read duals
  if (j) {
    VecReader<double> vr(f, binary, j);
    Handler().OnDualSolution(vr);
    if (!CheckReader( vr, readresult_ ))
      return readresult_;
  }

	i = nsv = have_options ? (int)z[3] : NumVars();
//  objno = ac_->no >= 0 ? 0 : -1;
  sstatus_seen = 0;
  if (binary) {
    if (!fread((char *)&L, sizeof(uiolen), 1, f) || L != L1)
      return ReportBadFormat();
    L1 = i * sizeof(real);
    if (!fread((char *)&L, sizeof(uiolen), 1, f) || L != L1)
      return ReportBadFormat();

    // Read primal solution
    if (i) {
      VecReader<double> vr(f, binary, i);
      Handler().OnPrimalSolution(vr);
      if (!CheckReader( vr, readresult_ ))
        return readresult_;
    }

    if (!fread((char *)&L, sizeof(uiolen), 1, f) || L != L1)
      return ReportBadFormat();

    /* do we have an objno ? */
    if (fread((char *)&L, sizeof(uiolen), 1, f)) {
      if (L != 2*sizeof(Long) && L != sizeof(integer))
        return ReportBadFormat();
      ui = (size_t) L;
      if (fread(Objno, ui, 1, f) != 1
       || fread(&L1, sizeof(uiolen), 1, f) != 1
       || L != L1)
        return ReportBadFormat();
      objno = (int)Objno[0];

      /* Submit objno and solve_code to Handler. */
      Handler().OnObjno(objno);
      Handler().OnSolveCode(Objno[1]);

      if (L == 2*sizeof(integer)) {
        switch(bsufread(f)) {
        case NLW2_SOLRead_Bad_Suffix:
badsuftable:
          serror("Bad suffix in '%s'\n", stub_);
          return NLW2_SOLRead_Bad_Suffix;
        case NLW2_SOLRead_OK:
          break;
        default:
          return readresult_;
        }
      }
      while(fread(&L, sizeof(uiolen), 1, f)) {}
    }
  }
  else {
    // Read primal solution
    if (i) {
      VecReader<double> vr(f, binary, i);
      Handler().OnPrimalSolution(vr);
      if (!CheckReader( vr, readresult_ ))
        return readresult_;
    }

    if (!fgets(buf,sizeof(buf), f))
      goto f_done;
    if (strncmp(buf,"objno ",6)) {
bad_objno:
      return ReportBadLine(buf);
      }
    x = strtod(s = buf+6, &se);
    if (se <= s)
      goto bad_objno;
    objno = (int)x;
    x = strtod(s = se, &se);
    if (se <= s)
      goto f_done;
    Objno[1] = (Long)x;

    /* Submit objno and solve_code to Handler. */
    Handler().OnObjno(objno);
    Handler().OnSolveCode(Objno[1]);

    switch(gsufread(f)) {
    case NLW2_SOLRead_Bad_Suffix:
      goto badsuftable;
    case NLW2_SOLRead_OK:
      break;
    default:
      return readresult_;
    }
  }
f_done:
  internal_rv_ = 0;
  return NLW2_SOLRead_OK;
}


/// Suffix head
struct SufHead {
  char sufid[8];
  int kind;
  int n;
  int namelen;
  int tablen;
};

/// Suffix read data
struct SufRead {
  SufHead h;
  char *name;
  char *tabname;
  char *table;
//  int *ip;
//  double *rp;
  std::vector<char> xp;
  int tablines;
//  int k;
//  int nmax;
};

template <class SOLHandler>
NLW2_SOLReadResultCode SOLReader2<SOLHandler>::bsufread(FILE* f) {
  uiolen L, L1;

  while(fread(&L, sizeof(uiolen), 1, f)) {
    SufRead SR;
    if (L < sizeof(SufHead))
      return NLW2_SOLRead_Bad_Suffix;
    if (fread(&SR.h, sizeof(SufHead), 1, f) != 1)
      return ReportEarlyEof();
    SR.tablines = SR.h.tablen - 1;
    if (strncmp(SR.h.sufid, "\nSuffix\n", 8)
        || sufheadcheck(&SR))
      return NLW2_SOLRead_Bad_Suffix;
    if (fread(SR.name, (size_t)SR.h.namelen, 1, f) != 1)
      return ReportEarlyEof();
    if (SR.h.tablen && fread(SR.table, (size_t)SR.h.tablen,
                             1, f) != 1)
      return NLW2_SOLRead_Bad_Suffix;
    SuffixInfo si(SR.h.kind, SR.name, SR.table);
    if (SR.h.kind & 4) {        // real-valued
      SuffixReader<double> sr(std::move(si), f, 1, SR.h.n);
      Handler().OnDblSuffix(sr);
      if (!CheckReader( sr, readresult_ ))
        return readresult_;
    } else {                    // int-valued
      SuffixReader<int> sr(std::move(si), f, 1, SR.h.n);
      Handler().OnIntSuffix(sr);
      if (!CheckReader( sr, readresult_ ))
        return readresult_;
    }
    if (!fread(&L1, sizeof(uiolen), 1, f) || L != L1)
      return ReportEarlyEof();
    // sufput(&SR, ac, newsufs);
  }
  return NLW2_SOLRead_OK;
}

/// Parse int
inline int
Lget(char **sp, int *Lp)
{
  int c;
  int L = 0;
  char *s = *sp;
  while(*s == ' ')
    s++;
  if ((c = *s++) < '0' || c > '9')
    return 1;
  L = c - '0';
  while((c = *s) >= '0' && c <= '9') {
    L = 10*L + c - '0';
    s++;
  }
  *Lp = L;
  *sp = s;
  switch(*s) {
  case '\r':
    if (*++s != '\n')
      break;
    *sp = s;
    /* no break */
  case 0:
  case ' ':
  case '\n':
    return 0;
  }
  return 1;
}


template <class SOLHandler>
NLW2_SOLReadResultCode SOLReader2<SOLHandler>::gsufread(FILE* f) {
  char *s, *se;
  size_t L;
  char buf[512];

  while(fgets(buf, sizeof(buf)-1, f)) {
    SufRead SR;
    if (strncmp(buf, "suffix ", 7))
      return ReportBadLine(buf);
    s = buf + 7;
    if (Lget(&s, &SR.h.kind)
        || Lget(&s, &SR.h.n)
        || Lget(&s, &SR.h.namelen)
        || Lget(&s, &SR.h.tablen)
        || Lget(&s, &SR.tablines))
      return ReportBadLine(buf);
    if (sufheadcheck(&SR))
      return ReportBadLine(buf);
    if (!fgets(buf, sizeof(buf)-1, f)
        || (buf[SR.h.namelen-1] != '\n'
            && (buf[SR.h.namelen-1] != '\r'
                || buf[SR.h.namelen] != '\n')))
      return ReportBadLine(buf);
    buf[SR.h.namelen-1] = 0;
    strcpy(SR.name, buf);
    if (SR.h.tablen) {
      s = SR.table;
      se = s + SR.h.tablen;
      for(int i = 1; i < SR.tablines; i++) {
        if (!fgets(s, se-s, f))
          return ReportEarlyEof();
        s += strlen(s);
      }
      if (!fgets(buf, sizeof(buf)-1, f))
        return ReportEarlyEof();
      if (!(L = strlen(buf)) || buf[--L] != '\n'
          || L >= (size_t)(se - s))
        return ReportBadLine(buf);
      if (L) {
        if (buf[L-1] != '\r' || --L) {
          buf[L] = 0;
          memcpy(s, buf, L);
        }
      }
    }
    SuffixInfo si(SR.h.kind, SR.name, SR.table);
    if (SR.h.kind & 4) {        // real-valued
      SuffixReader<double> sr(std::move(si), f, binary, SR.h.n);
      Handler().OnDblSuffix(sr);
      if (!CheckReader( sr, readresult_ ))
        return readresult_;
    } else {                    // int-valued
      SuffixReader<int> sr(std::move(si), f, binary, SR.h.n);
      Handler().OnIntSuffix(sr);
      if (!CheckReader( sr, readresult_ ))
        return readresult_;
    }
    //		sufput(&SR, ac, newsufs);
  }
  return NLW2_SOLRead_OK;
}

template <class SOLHandler>
int SOLReader2<SOLHandler>::sufheadcheck(SufRead* sr) {
  int n;

  n = (int)sr->h.n;
  if (sr->h.kind < 0 || sr->h.kind > 15 || n < 0 || sr->h.namelen < 2
   || sr->h.tablen < 0)
    return 1;
  i = (int)sr->h.kind & 3;
  if (sr->h.tablen
   && (sr->tablines > sr->h.tablen + 1 || sr->tablines < 1))
    return 1;
  sr->xp.resize((sr->h.tablen + 2*sr->h.namelen + 6));
  sr->name = (char*)sr->xp.data();
  sr->table = sr->name + sr->h.namelen;
  sr->tabname = sr->table + sr->h.tablen;
  return 0;
}

template <class SOLHandler>
NLW2_SOLReadResultCode SOLReader2<SOLHandler>::ReportEarlyEof() {
  serror("error reading '%s' (errno=%d)", stub_, errno);
  return readresult_ = NLW2_SOLRead_Early_EOF;
}

template <class SOLHandler>
NLW2_SOLReadResultCode SOLReader2<SOLHandler>::ReportBadFormat() {
  serror("Bad %s solution file '%s' (errno=%d)",
         bkind[binary], stub_, errno);
  return readresult_ = NLW2_SOLRead_Bad_Format;
}

template <class SOLHandler>
NLW2_SOLReadResultCode SOLReader2<SOLHandler>::
ReportBadLine(const std::string& line) {
  serror("Bad line in '%s': %s", stub_, line.c_str());
  return readresult_ = NLW2_SOLRead_Bad_Line;
}


template <class SOLHandler>
void SOLReader2<SOLHandler>::serror(
    const char* format, ...) {
  va_list args;
  va_start (args, format);
  auto n = std::vsnprintf (0, 0, format, args);
  va_end (args);
  va_start (args, format);
  std::vector<char> buf(n+1);
  auto n1 = std::vsnprintf (buf.data(), n+1, format, args);
  assert(n==n1);
  va_end (args);
  if (err_msg_.size())
    err_msg_ += '\n';
  err_msg_ += buf.data();
}



}  // namespace mp

#endif // SOLREADER2_HPP
