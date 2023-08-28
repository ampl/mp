/*
 A mathematical optimization solver.

 Copyright (C) 2012 AMPL Optimization Inc

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

#include <cctype>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>

#include <algorithm>
#include <stack>

#ifndef _WIN32
# include <strings.h>
# include <unistd.h>
# define MP_WRITE write
#else
#define NOMINMAX
# include <io.h>
# include <process.h>
# include <windows.h>
# define MP_WRITE _write
# define strcasecmp _stricmp
# undef max
#endif

#include "mp/common.h"
#include "mp/solver.h"

#include "mp/os.h"
#include "mp/utils-file.h"
#include "mp/utils-string.h"
#include "mp/clock.h"
#include "mp/rstparser.h"

#include "mp/ampls-cpp-api.h"

namespace {

const char *SkipSpaces(const char *s) {
  while (*s && isspace(*s))
    ++s;
  return s;
}

const char *SkipNonSpaces(const char *s) {
  while (*s && !isspace(*s))
    ++s;
  return s;
}

const char* SkipToEnd(const char* s) {
  while (*s && (*s != '\n'))
    ++s;
  return s;
}

const char* SkipToMatchingQuote(const char* s) {
  assert((*s == '\'') || (*s == '"'));
  char quote = s[0];
  ++s;
  while (*s != quote)
    ++s;
  return ++s;
}

struct Deleter {
  void operator()(mp::SolverOption* p) const { delete p; }
};

// A reStructuredText (RST) formatter.
class RSTFormatter : public rst::ContentHandler {
 private:
  fmt::Writer &writer_;
  mp::ValueArrayRef values_;
  std::stack<int> indents_;
  int indent_;
  int pos_in_line_;
  bool end_block_;  // true if there was no text after the last end of block

  enum {
    LITERAL_BLOCK_INDENT = 3,
    LIST_ITEM_INDENT     = 2
  };

  void Indent() {
    if (end_block_) {
      end_block_ = false;
      EndLine();
    }
    for (; pos_in_line_ < indent_; ++pos_in_line_)
      writer_ << ' ';
  }

  // Ends the current line.
  void EndLine() {
    writer_ << '\n';
    pos_in_line_ = 0;
  }

  // Writes a string doing a text reflow.
  void Write(fmt::CStringRef s);

 public:
  RSTFormatter(fmt::Writer &w, mp::ValueArrayRef values, int indent)
  : writer_(w), values_(values), indent_(indent),
    pos_in_line_(0), end_block_(false) {}

  void StartBlock(rst::BlockType type);
  void EndBlock();

  void HandleText(const char *text, std::size_t size);
  void HandleDirective(const char *type);
};

void RSTFormatter::Write(fmt::CStringRef s) {
  enum { MAX_LINE_LENGTH = 78 };
  const char *p = s.c_str();
  for (;;) {
    // Skip leading spaces.
    while (*p == ' ')
      ++p;
    const char *word_start = p;
    while (*p != ' ' && *p != '\n' && *p)
      ++p;
    const char *word_end = p;
    if (pos_in_line_ + (word_end - word_start) +
        (pos_in_line_ == 0 ? 0 : 1) > MAX_LINE_LENGTH) {
      EndLine();  // The word doesn't fit in the current line, start a new one.
    }
    if (pos_in_line_ == 0) {
      Indent();
    } else {
      // Separate words.
      writer_ << ' ';
      ++pos_in_line_;
    }
    std::size_t length = word_end - word_start;
    writer_ << fmt::StringRef(word_start, length);
    pos_in_line_ += static_cast<int>(length);
    if (*p == '\n') {
      EndLine();
      ++p;
    }
    if (!*p) break;
  }
}

void RSTFormatter::StartBlock(rst::BlockType type) {
  indents_.push(indent_);
  if (type == rst::LITERAL_BLOCK) {
    indent_ += LITERAL_BLOCK_INDENT;
  } else if (type == rst::LIST_ITEM) {
    Write("*");
    indent_ += LIST_ITEM_INDENT;
  }
}

void RSTFormatter::EndBlock() {
  indent_ = indents_.top();
  indents_.pop();
  end_block_ = true;
}

void RSTFormatter::HandleText(const char *text, std::size_t size) {
  std::string str(text, size);
  // Replace "``" with "\"".
  const char QUOTES[] = "``";
  std::size_t pos = 0;
  while ((pos = str.find(QUOTES, pos)) != std::string::npos) {
    str.replace(pos, sizeof(QUOTES) - 1, 1, '"');
    ++pos;
  }
  Write(str);
  size = writer_.size();
  if (size != 0 && writer_.data()[size - 1] != '\n')
    EndLine();
}

void RSTFormatter::HandleDirective(const char *type) {
  if (std::strcmp(type, "value-table") != 0)
    throw mp::Error("unknown directive {}", type);
  if (values_.size() == 0)
    throw mp::Error("no values to format");
  std::size_t max_len = 0;
  for (mp::ValueArrayRef::iterator
      i = values_.begin(), e = values_.end(); i != e; ++i) {
    max_len = std::max(max_len, std::strlen(i->value));
  }

  // If the values don't have descriptions indent them as list items.
  if (!values_.begin()->description)
    indent_ += LIST_ITEM_INDENT;
  for (mp::ValueArrayRef::iterator
      i = values_.begin(), e = values_.end(); i != e; ++i) {
    Indent();
    writer_ << fmt::pad(i->value, static_cast<unsigned>(max_len));
    if (i->description) {
      static const char SEP[] = " -";
      writer_ << SEP;
      int saved_indent = indent_;
      indent_ += static_cast<int>(max_len + sizeof(SEP));
      pos_in_line_ = indent_;
      Write(i->description);
      indent_ = saved_indent;
    }
    EndLine();
  }

  // Unindent values if necessary and end the block.
  if (!values_.begin()->description)
    indent_ -= LIST_ITEM_INDENT;
  end_block_ = true;
}

}  // namespace

namespace mp {

namespace internal {

void FormatRST(fmt::Writer &w,
    fmt::CStringRef s, int indent, ValueArrayRef values) {
  RSTFormatter formatter(w, values, indent);
  rst::Parser parser(&formatter);
  parser.Parse(s.c_str());
}

int OptionHelper<int>::Parse(const char *&s, bool) {
  char *end = 0;
  long value = std::strtol(s, &end, 10);
  s = end;
  return value;
}

double OptionHelper<double>::Parse(const char *&s, bool) {
  char *end = 0;
  double value = std::strtod(s, &end);
  s = end;
  return value;
}
bool quoted(const char* s) {
  return (*s == '\'' || *s == '"');
}
std::string OptionHelper<std::string>::Parse(const char *&s, bool splitString) {
  const char *start = s;
  if (splitString) // if the string has been already split (by the command line parser)
  {
    s = SkipToEnd(s);
    return std::string(start, s - start);
  }
  if (quoted(s))
  {
    s = SkipToMatchingQuote(s);
    return std::string(start + 1, s - start - 2);
  }
  else
  {
    s = SkipNonSpaces(s);
    return std::string(start, s - start);
  }

}

SolverAppOptionParser::SolverAppOptionParser(BasicSolver &s)
  : solver_(s), echo_solver_options_(true) {
  // Add standard command-line options.
  OptionList::Builder<SolverAppOptionParser> app_options(options_, *this);
  app_options.Add<&SolverAppOptionParser::ShowUsage>(
        '?', "show usage and exit");
  app_options.Add<&SolverAppOptionParser::EndOptions>('-', "end of options");
  app_options.AddWithParam<&SolverAppOptionParser::ShowSolverOptions>(
        '=', "show solver options and exit");
  app_options.Add<&SolverAppOptionParser::ShowSolverOptionsASL>(
    'a', "show solver options (ASL style, 1st synonyms if provided) and exit");
  app_options.Add<&SolverAppOptionParser::DontEchoSolverOptions>(
        'e', "suppress echoing of assignments");
  app_options.Add<&SolverAppOptionParser::WantSol>(
        's', "write .sol file (without -AMPL)");
  OptionList::Builder<mp::BasicSolver> options(options_, s);
  options.Add<&mp::BasicSolver::ShowVersion>('v', "show version and exit");
  options.Add<&mp::BasicSolver::ShowConstraintDescriptions>('c',
                                             "show constraint descriptions and exit");
  // TODO: if solver supports functions add options -ix and -u
  // "ix", "import user-defined functions from x; -i? gives details",
  // "u",  "just show available user-defined functions"
}

bool SolverAppOptionParser::ShowUsage() {
  solver_.Print("usage: {} [options] stub [-AMPL] [<assignment> ...]\n",
                solver_.name());
  solver_.Print("\nOptions:\n");
  for (OptionList::iterator
       i = options_.begin(), end = options_.end(); i != end; ++i) {
    solver_.Print("\t-{}  {}\n", i->name, i->description);
  }
  return false;
}


bool SolverAppOptionParser::ShowSolverOptionsASL() {
  struct OptionASLNameLess {
    bool operator()(
      const SolverOption* lhs, const SolverOption* rhs) const {
      int cmp = strcasecmp(lhs->name_ASL(), rhs->name_ASL());
      return cmp < 0;
    }
  };
  fmt::MemoryWriter writer;
  const char* option_header = solver_.option_header();
  internal::FormatRST(writer, option_header);
  if (*option_header)
    writer << '\n';
  solver_.Print("{}", writer.c_str());
  solver_.Print("Options:\n");
  const int DESC_INDENT = 6;

  // Create set with ASL ordering
  std::set<const SolverOption*, OptionASLNameLess> set;
  for (Solver::option_iterator
    i = solver_.option_begin(), end = solver_.option_end(); i != end; ++i) {
    auto cc = (&(*i));
    set.insert(cc);
  }
  for (std::set<const SolverOption*, OptionASLNameLess>::const_iterator i = set.begin();
    i != set.end(); ++i)
  {
    writer.clear();
    writer << '\n' << (*i)->name_ASL() << '\n';
    (*i)->format_description(writer, DESC_INDENT);
    solver_.Print("{}", fmt::StringRef(writer.data(), writer.size()));
  }
  return false;
}
bool contains(const char* name, const char* substr) {
  std::string n(name);
  std::string sub = fmt::format("{}:", substr);
  return n.find(sub) != std::string::npos;


}
bool SolverAppOptionParser::ShowSolverOptions(const char* param) {
  fmt::MemoryWriter writer;
  const char* option_header = solver_.option_header();
  internal::FormatRST(writer, option_header);
  if (*option_header)
    writer << '\n';
  solver_.Print("{}", writer.c_str());
  if(param)
    solver_.Print("{} ", param);
  solver_.Print("Options:\n");
  const int DESC_INDENT = 6;
  for (Solver::option_iterator
    i = solver_.option_begin(), end = solver_.option_end(); i != end; ++i) {
    if (param && strlen(param)>0)
      if (!contains(i->name(), param))
        continue;
    writer.clear();
    writer << '\n' << i->name();
    const auto& syns = i->inline_synonyms();
    if (auto sz = syns.size()) {
      writer << " (";
      for (const auto& s : syns) {
        writer << s;
        if (--sz)
          writer << ", ";
      }
      writer << ')';
    }
    writer << '\n';
    i->format_description(writer, DESC_INDENT);
    solver_.Print("{}", fmt::StringRef(writer.data(), writer.size()));
  }
  return false;
}

const char *SolverAppOptionParser::Parse(char **&argv) {
  ++argv;
  char opt = ParseOptions(argv, options_);
  if (opt && opt != '-') return 0;
  const char *stub = *argv;
  if (!stub) {
    ShowUsage();
    return 0;
  }
  ++argv;
  if (*argv && std::strcmp(*argv, "-AMPL") == 0) {
    solver_.set_ampl_flag();
    solver_.set_wantsol(1);
    ++argv;
  }
  return stub;
}

mp::internal::atomic<const char*> SignalHandler::signal_message_ptr_;
mp::internal::atomic<unsigned> SignalHandler::signal_message_size_;
mp::internal::atomic<InterruptHandler> SignalHandler::handler_;
mp::internal::atomic<void*> SignalHandler::data_;

volatile std::sig_atomic_t SignalHandler::stop_ = 1;

#ifdef _WIN32
// Signal repeater for Windows.
// It is required since GenerateConsoleCtrlEvent doesn't
// work reliably with CTRL_C_EVENT across processes.
void RunSignalRepeater(void *arg) {
  using mp::internal::SignalRepeater;
  SignalRepeater *repeater = reinterpret_cast<SignalRepeater*>(arg);
  HANDLE in = reinterpret_cast<HANDLE>(repeater->in());
  HANDLE out = reinterpret_cast<HANDLE>(repeater->out());
  for (;;) {
    struct {
      int sig;
      int pid;
      int num;
    } data = {0, 0, 0};
    DWORD count = 0;
    if (!ReadFile(in, &data, sizeof(data), &count, 0) || count == 0)
      break;
    if (data.pid > 0 && data.pid != GetCurrentProcessId())
      data.sig = 0;
    int sig = 0;
    // Translate signal codes.
    enum { SigInt = 2, SigTerm = 15, SigBreak = 21 };
    switch (data.sig) {
    case SigInt:
      sig = SIGINT;
      break;
    case SigBreak:
      sig = SIGBREAK;
      break;
    case SigTerm:
      sig = SIGTERM;
      break;
    }
    if (sig != 0) {
      data.num++;
      if (data.pid >= 0)
        data.sig = 0;
    }
    WriteFile(out, &data, count, &count, 0);
    if (sig)
      std::raise(sig);
    Sleep(50);
  }
}

mp::internal::SignalRepeater::SignalRepeater(const char *s) : in_(0), out_(0) {
  if (!s)
    return;
  char *end = 0;
  in_ = _strtoui64(s, &end, 10);
  if (in_ == 0 || *end != ',')
    return;
  s = end + 1;
  out_ = _strtoui64(s, &end, 10);
  if (out_ != 0)
    _beginthread(RunSignalRepeater, 0, this);
}
#endif

SignalHandler::SignalHandler(BasicSolver &s)
  : solver_(s), message_(fmt::format("\n<BREAK> ({})\n", s.name())),
    repeater_(std::getenv("SW_sigpipe")) {
  solver_.set_interrupter(this);
  signal_message_ptr_ = message_.c_str();
  signal_message_size_ = static_cast<unsigned>(message_.size());
  std::signal(SIGINT, HandleSigInt);
  stop_ = 0;
}

SignalHandler::~SignalHandler() {
  solver_.set_interrupter(0);
  stop_ = 1;
  handler_ = 0;
  signal_message_size_ = 0;
}

void SignalHandler::SetHandler(InterruptHandler handler, void *data) {
  handler_ = handler;
  data_ = data;
}

void SignalHandler::HandleSigInt(int sig) {
  unsigned count = 0;
  do {
    // Use asynchronous-safe function write instead of printf!
    int result = MP_WRITE(1, signal_message_ptr_ + count,
        signal_message_size_ - count);
    if (result < 0) break;
    count += result;
  } while (count < signal_message_size_);
  if (stop_) {
    // Use asynchronous-safe function _exit instead of exit!
    _exit(1);
  }
  stop_ = 1;
  if (InterruptHandler handler = handler_)
    handler(data_);
  // Restore the handler since it might have been reset before the handler
  // is called (this is implementation defined).
  std::signal(sig, HandleSigInt);
}


void PrintSolution(const double *values, int num_values, const char *name_col,
                   const char *value_col, NameProvider &np) {
  if (!values || num_values == 0)
    return;
  std::size_t name_len = std::strlen(name_col);
  std::size_t name_field_width = name_len;
  for (int i = 0; i < num_values; ++i)
    name_field_width = std::max(name_field_width, np.name(i).size());
  name_field_width += 2;
  fmt::printf("\n%-*s%s\n", name_field_width, name_col, value_col);
  for (int i = 0; i < num_values; ++i) {
    double value = values[i];
    fmt::printf("%-*s%.17g\n", name_field_width, np.name(i), value ? value : 0);
  }
}
}  // namespace internal

bool Solver::OptionNameLess::operator()(
    const SolverOption *lhs, const SolverOption *rhs) const {
  int cmp = strcasecmp(lhs->name(), rhs->name());
  return cmp < 0;
}


BasicSolver::BasicSolver() :
  name_("dummy"), long_name_("long_dummy") { }

BasicSolver::BasicSolver(
    fmt::CStringRef name, fmt::CStringRef long_name, long date, int flags) {
  InitMetaInfoAndOptions(name, long_name, date, flags);
}

void BasicSolver::InitMetaInfoAndOptions(
    fmt::CStringRef name, fmt::CStringRef long_name, long date, int flags) {
  name_ = name.c_str();
  long_name_ = (long_name.c_str() ? long_name : name).c_str();
  date_ = date;

  version_ = long_name_;

  struct VersionOption : SolverOption {
    BasicSolver &s;
    VersionOption(BasicSolver &s) : SolverOption("tech:version version",
        "Single-word phrase: report version details "
        "before solving the problem.", ValueArrayRef(), true), s(s) {}

    void Write(fmt::Writer &w) { w << ((s.bool_options_ & SHOW_VERSION) != 0); }
    void Parse(const char *&, bool) { s.bool_options_ |= SHOW_VERSION; }
    Option_Type type() {
      return Option_Type::BOOL;
    }
  };
  AddOption(OptionPtr(new VersionOption(*this)));

  AddStrOption(
        "tech:optionfile optionfile option:file",
        "Name of solver option file. "
        " (surrounded by 'single' or "
        "\"double\" quotes if the name contains blanks). "
        "Lines that start with # are ignored.  Otherwise, each nonempty "
        "line should contain \"name=value\".",
        &Solver::GetOptionFile, &Solver::UseOptionFile);


  AddIntOption(
        "tech:wantsol wantsol",
        "In a stand-alone invocation (no ``-AMPL`` on the command line), "
        "what solution information to write.  Sum of\n"
        "\n"
        "| 1 - Write ``.sol`` file\n"
        "| 2 - Primal variables to stdout\n"
        "| 4 - Dual variables to stdout\n"
        "| 8 - Suppress solution message.",
        &Solver::GetWantSol, &Solver::SetWantSol);

  AddIntOption(
        "obj:no objno",
        "Objective to optimize:\n"
        "\n"
        "| 0 - None\n"
        "| 1 - First (default, if available)\n"
        "| 2 - Second (if available), etc.\n",
        &Solver::GetObjNo, &Solver::SetObjNo);

  struct BoolOption : TypedSolverOption<int> {
    bool &value_;
    BoolOption(bool &value, const char *name, const char *description)
    : TypedSolverOption<int>(name, description), value_(value) {}

    void GetValue(fmt::LongLong &value) const { value = value_; }
    void SetValue(fmt::LongLong value) {
      if (value != 0 && value != 1)
        throw InvalidOptionValue(name(), value);
      value_ = value != 0;
    }
  };


  AddOption(OptionPtr(new BoolOption(debug_, "tech:debug debug",
                                     "0*/1: whether to assist testing & debugging, e.g., "
                                     "by outputting auxiliary information.")));

  if ((flags & MULTIPLE_OBJ) != 0) {
    AddOption(OptionPtr(new BoolOption(multiobj_, "obj:multi multiobj",
                                       "0*/1:  Whether to use multi-objective optimization. "
                                       "If set to 1 multi-objective optimization is performed using "
                                       "lexicographic method with the first objective treated as the most "
                                       "important, then the second objective and so on.")));
  }

  AddOption(OptionPtr(new BoolOption(timing_, "tech:timing timing",
      "0*/1: Whether to display timings for the run.")));

  if ((flags & MULTIPLE_SOL) != 0) {

    AddOption(OptionPtr(new BoolOption(count_solutions_, "sol:count countsolutions",
        "0*/1: Whether to count the number of solutions "
        "and return it in the ``.nsol`` problem suffix.")));

    AddStrOption(
          "sol:stub solstub solutionstub",
          "Stub for solution files.  If ``solutionstub`` is specified, "
          "found solutions are written to files (``solutionstub & '1' & "
          "'.sol'``) ... (``solutionstub & Current.nsol & '.sol'``), where "
          "``Current.nsol`` holds the number of returned solutions.  That is, "
          "file names are obtained by appending 1, 2, ... ``Current.nsol`` to "
          "``solutionstub``.",
          &Solver::GetSolutionStub, &Solver::SetSolutionStub);
  }
}

/// Process lines with a custom line processor.
/// Used to parse an options file.
static void ProcessLines_AvoidComments(std::istream& stream,
                  std::function<void(const char*)> processor) {
  std::string line;
  while (stream.good() && !stream.eof()) {
    std::getline(stream, line);
    if (line.size()) {
      auto itfirstns = std::find_if(line.begin(), line.end(),
                                  [](char c){ return !std::isspace(c); });
      if (line.end()!=itfirstns &&
              '#'!=*itfirstns) {                  // Skip commented line
        processor(line.c_str()+(itfirstns-line.begin()));
      }
    }
  }
}

void BasicSolver::UseOptionFile(const SolverOption &, fmt::StringRef value) {
  option_file_save_ = value;
  std::ifstream ifs(value);
  if (ifs.good())
    ProcessLines_AvoidComments(ifs,
                               [this](const char* s)
    { ParseOptionString(s, option_flag_save_); });
  if (!ifs.good() && !ifs.eof())
    MP_RAISE(fmt::format("Failed to read option file '{}': {}",
                         value, std::strerror(errno)));
}

SolverOptionManager::~SolverOptionManager() {
  std::for_each(options_.begin(), options_.end(), Deleter());
}

void BasicSolver::AddWarning(
    std::string key, std::string msg, bool replace) {
  auto& v = GetWarningsMap()[ std::move(key) ];
  if (!v.first++      // only remember the 1st detailed message
      || replace)     // unless asked to replace
    v.second = std::move(msg);
}

std::string BasicSolver::GetWarnings() const {
  if (GetWarningsMap().size()) {
    std::string wrn = "------------ WARNINGS ------------\n";
    for (const auto& e: GetWarningsMap())
      wrn += ToString(e) + '\n';
    return wrn;
  }
  return "";
}

void BasicSolver::PrintWarnings() {
  auto warnings = GetWarnings();
  if (warnings.size())
    Print('\n' + warnings);
}

std::string BasicSolver::ToString(
    const WarningsMap::value_type& wrn) {
  return fmt::format(
        "WARNING.  {} case(s) of \"{}\". One of them:\n  {}",
        wrn.second.first, wrn.first, wrn.second.second);
}

#ifdef MP_DATE
bool BasicSolver::ShowVersion() {
  Print("{} ({})", version_, MP_SYSINFO);
  if (date_ > 0)
    Print(", driver({})", date_);
  Print(", MP({})\n", MP_DATE);
  if (!license_info_.empty())
    Print("{}\n", license_info_);
  if (!this->set_external_libs().empty())
    Print("External libraries:\n{}", this->set_external_libs());
  return false;
}
#else
// For use with ASL.
extern "C" char sysdetails_ASL[];
extern "C" const char *Lic_info_ASL, *Version_Qualifier_ASL;
extern "C" long ASLdate_ASL;
extern "C" void Mach_ASL();
bool Solver::ShowVersion() {
  if (*Version_Qualifier_ASL) {
    Mach_ASL(); // may adjust Version_Qualifier_ASL
    Print("{}", Version_Qualifier_ASL);
  }
  const char *sysdetails = sysdetails_ASL;
  Print("{} ({})", version_, sysdetails);
  if (date_ > 0)
    Print(", driver({})", date_);
  Print(", ASL({})\n", ASLdate_ASL);
  if (Lic_info_ASL && *Lic_info_ASL)
    Print("{}\n", Lic_info_ASL);
  return false;
}
#endif

bool BasicSolver::ShowConstraintDescriptions() {
  Print("{}\n\n", constr_descr_header_);
  if (constr_descr_.empty())
    Print("No constraint descrptions filled.\n");
  else {
    int i=0;
    for (const auto& cd: constr_descr_) {
      Print("{:3}.  {:20}:  {}\n", ++i, cd.first, cd.second);
    }
  }
  return false;
}

SolverOption::SolverOption(
    const char *names_list, const char *description,
    ValueArrayRef values, bool is_flag)
: description_(description),
  values_(values), is_flag_(is_flag) {
  auto synonyms = split_string(names_list);
  if (synonyms.empty())
    throw std::logic_error("Empty option name list");
  name_ = synonyms.front();
  for (size_t i=1; i<synonyms.size(); ++i)
    inline_synonyms_.push_back(synonyms[i]);
  /// Wildcards
  auto wc_pos = name_.find_first_of('*');
  if (std::string::npos != wc_pos) {
    wc_headtails_.push_back(wc_split(name_));
    for (const auto& syn: inline_synonyms_)
      wc_headtails_.push_back(wc_split(syn));
  }
}

SolverOption::WCHeadTail SolverOption::wc_split(
    const std::string &name) {
  assert(name.size()>1);
  auto wc_pos = name.find_first_of('*');
  assert((std::string::npos != wc_pos));
  assert(wc_pos == name.find_last_of('*'));
  return {
    name.substr(0, wc_pos),
    name.substr(wc_pos+1, std::string::npos)
  };
}

bool SolverOption::wc_match(const std::string &key) {
  for (const auto& wcht: wc_headtails_) {
    if (0==key.rfind(wcht.first, 0) &&
        key.size()>wcht.second.size() &&
        key.size()-wcht.second.size() ==
          key.rfind((wcht.second))) {
      wc_key_last_ = key;
      wc_body_last_ = key.substr(
            wcht.first.size(), key.size()-wcht.second.size()-wcht.first.size());
      return true;
    }
  }
  return false;
}

SolverOption *SolverOptionManager::FindOption(
    const char *name, bool wildcardvalues) const {
  struct DummyOption : SolverOption {
    DummyOption(const char *name) : SolverOption(name, "") {}
    void Write(fmt::Writer &) {}
    void Parse(const char *&, bool) {}
    Option_Type type() {
      return Option_Type::BOOL;
    }
  };
  DummyOption option(name);
  // find by name
  OptionSet::const_iterator i = options_.find(&option);
  if (i != options_.end()) {
    if ((*i)->is_wildcard() && wildcardvalues)
      return 0;
    return *i;
  }

  // find by inline synonyms and wildcards. Case-insensitive
  std::string name_str {name};
  for (OptionSet::const_iterator i = options_.begin();
    i != options_.end(); ++i)
  {
    if (std::find_if((*i)->inline_synonyms().begin(),
                    (*i)->inline_synonyms().end(),
                    [&name_str](const std::string& syn) {
                     return 0==strcasecmp(name_str.c_str(), syn.c_str()); } ) !=
        (*i)->inline_synonyms().end()) {
      return *i;
    }
    /// Wildcards
    if (wildcardvalues && (*i)->wc_match(name))
      return *i;
  }
  return 0;
}

void BasicSolver::ParseOptionString(
    const char *s, unsigned flags) {
  for (;;) {
    if (!*(s = SkipSpaces(s)))
      return;

    // Parse the option name.
    const char *name_start = s;
    while (*s && !std::isspace(*s) && *s != '=')
      ++s;
    fmt::internal::MemoryBuffer<char, 50> name;
    std::size_t name_size = s - name_start;
    name.resize(name_size + 1);
    for (std::size_t i = 0; i < name_size; ++i)
      name[i] = name_start[i];
    name[name_size] = 0;

    // Check if we have an '=' sign,
    // skip all spaces intil next token.
    bool equal_sign = false;
    s = SkipSpaces(s);
    if (*s == '=') {
      s = SkipSpaces(s + 1);
      equal_sign = true;
    }

    // Parse option name.
    SolverOption *opt = FindOption(&name[0], true);
    if (!opt) {
      HandleUnknownOption(&name[0]);
      continue;       // in case it does not throw
    }

    // If user asks the default/current value.
    if (*s == '?') {
      char next = s[1];
      if (!next || std::isspace(next)) {
        ++s;
        if ((flags & NO_OPTION_ECHO) == 0) {
          Print("  {}", opt->echo_with_value() + '\n');
        }
        continue;
      }
    }

    /// Parse value if needed.
    if (equal_sign) {
      // No '=' for flags.
      if (opt->is_flag()) {
        ReportError(
              "Option \"{}\" doesn't accept an argument",
              &name[0]);
        s = SkipNonSpaces(s);    // In case we go on
        continue;
      } else {
        opt->Parse(s, flags & FROM_COMMAND_LINE);
      }
    } else {
      if (opt->is_flag()) {         // Might set some flag
        opt->Parse(s, flags & FROM_COMMAND_LINE);
      } else {                      // Even w/o '=' sign
        opt->Parse(s, flags & FROM_COMMAND_LINE);
      }
    }

    // Echo name [= value].
    if ((flags & NO_OPTION_ECHO) == 0) {
      Print("  {}", opt->echo_with_value() + '\n');
    }
  }
}

bool BasicSolver::ParseOptions(char **argv, unsigned flags, const ASLProblem *) {
  has_errors_ = false;
  bool_options_ &= ~SHOW_VERSION;
  option_flag_save_ = flags;
  bool had_exe_name_option_var = false;
  /// Look for a <solver>_options env var.
  /// First try the executable name.
  const char *s = exe_path();
  if (std::strlen(s)) {
    path p(s);
    auto exe_basename = p.filename().string();
    auto pt = exe_basename.rfind('.');
    if (std::string::npos != pt) {
      auto ext = exe_basename.substr(pt);
      if (".exe"==ext || ".app"==ext)
        exe_basename = exe_basename.substr(0, pt);
    }
    if (const char *s = std::getenv((exe_basename + "_options").c_str())) {
      ParseOptionString(s, flags);
      had_exe_name_option_var = true;
    }
  }
  // Otherwise try '<solver_name>_options'
  if (!had_exe_name_option_var)
  if (const char *s = std::getenv((name_ + "_options").c_str())) {
    ParseOptionString(s, flags);
  }
  flags |= FROM_COMMAND_LINE;        // proceed to command-line options
  if (argv) {
    while (const char *s = *argv++) {
      ParseOptionString(s, flags);
    }
  }
  if ((bool_options_ & SHOW_VERSION) != 0)
    ShowVersion();
  return !has_errors_;
}

BasicSolver::DoubleFormatter BasicSolver::FormatObjValue(double value) {
  if (obj_precision_ < 0) {
    const char *s =  std::getenv("objective_precision");
    obj_precision_ = s ? std::atoi(s) : 0;
    if (obj_precision_ == 0)
      obj_precision_ = DEFAULT_PRECISION;
  }
  DoubleFormatter formatter = {value, obj_precision_};
  return formatter;
}

const char* SolverOption::name_ASL() const {
  return inline_synonyms_.size() == 0 ? name() :
    inline_synonyms_[0].c_str();
}
void SolverOption::add_synonyms_front(const char* names_list) {
  auto synonyms = split_string(names_list);
  inline_synonyms_.insert(inline_synonyms_.begin(),
                          synonyms.begin(), synonyms.end());
}
void SolverOption::add_synonyms_back(const char* names_list) {
  auto synonyms = split_string(names_list);
  inline_synonyms_.insert(inline_synonyms_.end(),
                          synonyms.begin(), synonyms.end());
}

void SolverOptionManager::AddOption(OptionPtr opt) {
  // First insert the option, then release a pointer to it. Doing the other
  // way around may lead to a memory leak if insertion throws an exception.
  if (!options_.insert(opt.get()).second)
    throw std::logic_error(
        fmt::format("Option {} already defined", opt.get()->name()));
  opt.release();
}

void SolverOptionManager::AddOptionSynonyms_Inline_Front(
    const char* names_list, const char* realName)
{
  SolverOption* real = FindOption(realName);
  if (!real)
    throw std::logic_error(
        fmt::format("Option {} referred to by synonyms {} is unknown",
                    realName, names_list));
  real->add_synonyms_front(names_list);
}

void SolverOptionManager::AddOptionSynonyms_Inline_Back(
    const char* names_list, const char* realName)
{
  SolverOption* real = FindOption(realName);
  if (!real)
    throw std::logic_error(
        fmt::format("Option {} referred to by synonyms {} is unknown",
                    realName, names_list));
  real->add_synonyms_back(names_list);
}

/// An "out-of-line" option synonym
class SolverOptionSynonym : public SolverOption
{
  SolverOption* real_;
  std::string desc_;
public:
  SolverOptionSynonym(const char* names, SolverOption& real) :
    SolverOption(names, ""), real_(&real) {
    desc_ = fmt::sprintf("Synonym for %s.", real_->name());
    set_description( desc_.c_str() );
  }
  SolverOption* getRealOption() const { return real_; }

  virtual void Write(fmt::Writer& w) {
    real_->Write(w);
  }
  virtual void Parse(const char*& s, bool b) {
    real_->Parse(s, b);
  }

  virtual std::string echo() {
    return fmt::format("{} ({})", name(), real_->echo());
  }
  Option_Type type() {
    return real_->type();
  }
};

void SolverOptionManager::AddOptionSynonyms_OutOfLine(
    const char* name, const char* realName)
{
  SolverOption* real = FindOption(realName);
  if (!real)
    throw std::logic_error(
        fmt::format("Option {} referred to by {} is unknown", realName, name));
  OptionPtr opt = OptionPtr(new SolverOptionSynonym(name, *real));
  options_.insert(opt.get());
  opt.release();
}


}  // namespace mp


//// AMPLS C/C++ API

/// Option that owns the formatted description
struct AMPLSOption {
  std::string name;
  std::string description;
  mp::SolverOption::Option_Type type;
  AMPLSOption(const char* name, const char* descr, mp::SolverOption::Option_Type type) {
    this->name = name;
    this->description = descr;
    this->type = type;
  }
};
/// Internal AMPLS API Infos
typedef struct AMPLS_MP__internal_T {
  /// The Backend
  std::unique_ptr<mp::BasicBackend> p_be_;
  /// An output handler. Can be a parameter
  mp::OutputHandler output_h_;
  /// Messages from Env::GetWarnings()
  std::vector<std::string> msg_wrn_;
  /// Messages in addition to Env::GetWarnings()
  std::vector<std::string> msg_extra_;
  /// The 0-terminated char*'s to the messages,
  /// filled by AMPLSGetMessages()
  std::vector<const char*> msg_pchar_;
  /// Array of options
  std::vector<AMPLSOption> options_;
  std::vector<AMPLS_C_Option> options_c_;
  std::string optionvalue_;
} AMPLS_MP__internal;


AMPLS_MP_Solver* AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend> p_be,
  CCallbacks cb = {}) {
  AMPLS_MP_Solver* slv = new AMPLS_MP_Solver();
  AMPLS__internal__TryCatchWrapper( slv, [&]() {
    auto ii = new AMPLS_MP__internal();
    slv->internal_info_ = ii;
    ii->p_be_ = std::move(p_be);
    auto be = ii->p_be_.get();
    be->GetCallbacks() = cb;
    be->set_output_handler(&ii->output_h_);
    char* argv[] = {(char*)"ampls-driver", nullptr};
    be->Init(argv);
    be->set_wantsol(1);       // user can still modify by 'wantsol=...'
    return 0;
  } );
  return slv;
}

void AMPLS__internal__Close(AMPLS_MP_Solver* slv) {
  assert(slv->internal_info_);
  delete (AMPLS_MP__internal*)slv->internal_info_;
  delete slv;
}

void AMPLSReadExtras(AMPLS_MP_Solver* slv) {
  auto be = AMPLSGetBackend(slv);
  be->InputExtras();
}

void AMPLSSolve(AMPLS_MP_Solver* slv) {
  auto be = AMPLSGetBackend(slv);
  be->Solve();
}
mp::BasicBackend* AMPLSGetBackend(AMPLS_MP_Solver* slv) {
  assert(slv->internal_info_);
  return ((AMPLS_MP__internal*)(slv->internal_info_))->p_be_.get();
}

/// Load model incl. suffixes
int AMPLSLoadNLModel(AMPLS_MP_Solver* slv,
                     const char* nl_filename,
                     char** options) {
  return AMPLS__internal__TryCatchWrapper( slv, [=]() {
    std::string nl_filename_ = nl_filename;
    auto filename_no_ext_ = nl_filename_;
    const char *ext = std::strrchr(nl_filename, '.');
    if (!ext || std::strcmp(ext, ".nl") != 0)
      nl_filename_ += ".nl";
    else
      filename_no_ext_.resize(filename_no_ext_.size() - 3);

    auto be = AMPLSGetBackend(slv);
    be->ReadNL(nl_filename, filename_no_ext_, options);
    be->InputExtras();
  } );
}

int AMPLSReportResults(AMPLS_MP_Solver* slv, const char* solFileName = NULL) {
  return AMPLS__internal__TryCatchWrapper( slv, [=]() {
    auto be = AMPLSGetBackend(slv);
    // If solFileName is NULL, we scrap any previous override
    be->OverrideSolutionFile(solFileName == NULL ? "" : solFileName);
    be->ReportResults();
  } );
}

void AMPLSAddMessage(AMPLS_MP_Solver* slv, const char* msg) {
  assert(slv->internal_info_);
  assert(msg);
  ((AMPLS_MP__internal*)(slv->internal_info_))->msg_extra_.
      push_back(msg);
}

const char* const * AMPLSGetMessages(AMPLS_MP_Solver* slv) {
  auto be = AMPLSGetBackend(slv);
  auto ii = ((AMPLS_MP__internal*)(slv->internal_info_));
  auto& msg_wrn = ii->msg_wrn_;
  msg_wrn.clear();
  const auto& msg_extra = ii->msg_extra_;
  auto& pchar_vec = ii->msg_pchar_;
  pchar_vec.clear();

  /// Add Backend warnings
  for (const auto& wrn: be->GetWarningsMap())
    msg_wrn.push_back(mp::BasicSolver::ToString(wrn) + '\n');

  /// Add Backend warnings
  for (const auto& wrn: msg_wrn)
    pchar_vec.push_back(wrn.c_str());

  /// Add extra messages
  for (const auto& msg: msg_extra)
    pchar_vec.push_back(msg.c_str());

  pchar_vec.push_back(nullptr);                // final 0
  return pchar_vec.data();
}


AMPLS_C_Option*  AMPLSGetOptions(AMPLS_MP_Solver* slv) {
  auto be = AMPLSGetBackend(slv);
  auto ii = ((AMPLS_MP__internal*)(slv->internal_info_));
  if (ii->options_.size() == 0)
  {
    auto end = be->option_end();
    for (auto it = be->option_begin(); it != end; ++it) {
      AMPLSOption o(it->name(), it->format_description(4).c_str(), it->type());
      ii->options_.push_back(o);
    }
    for (const AMPLSOption& o : ii->options_) {
      ii->options_c_.push_back({ o.name.c_str(), o.description.c_str(), (int)o.type});
    }
  }
  ii->options_c_.push_back({});
  return ii->options_c_.data();
}

template <typename T> int AMPLSetOption(AMPLS_MP_Solver* slv,
  const char* name, T v) {
  auto be = AMPLSGetBackend(slv);
  try {
    auto opt = be->GetOption(name);
    opt->SetValue(v);
    return 0;
  }
  catch (const mp::OptionError& o) {
    ((AMPLS_MP__internal*)(slv->internal_info_))->msg_extra_.
      push_back(o.what());
    return 1;
  }
  catch (const std::runtime_error& e) {
    ((AMPLS_MP__internal*)(slv->internal_info_))->msg_extra_.
      push_back(e.what());
  }
  catch (const std::exception& ) {
  }
  return 2;
}
template <typename T> int AMPLSGetOption(AMPLS_MP_Solver* slv,
  const char* name, T v) {
  auto be = AMPLSGetBackend(slv);
  try {
    auto opt = be->GetOption(name);
    opt->GetValue(*v);
    return 0;
  }
  catch (const mp::OptionError& o) {
    ((AMPLS_MP__internal*)(slv->internal_info_))->msg_extra_.
      push_back(o.what());
    return 1;
  }
  catch (const std::runtime_error& e) {
    ((AMPLS_MP__internal*)(slv->internal_info_))->msg_extra_.
      push_back(e.what());
  }
  catch (const std::exception&) {
  }
  return 2;
}
int AMPLSSetIntOption(AMPLS_MP_Solver* slv,
  const char* name, int v) {
  return AMPLSetOption(slv, name, v);
}
int AMPLSSetDblOption(AMPLS_MP_Solver* slv,
  const char* name, double v) {
  return AMPLSetOption(slv, name, v);
}
int AMPLSSetStrOption(AMPLS_MP_Solver* slv,
  const char* name, const char* v) {
  return AMPLSetOption(slv, name, v);
}

int AMPLSGetIntOption(AMPLS_MP_Solver* slv,
  const char* name, int* v) {
  return AMPLSGetOption(slv, name, v);
}
int AMPLSGetDblOption(AMPLS_MP_Solver* slv,
  const char* name, double* v) {
  return AMPLSGetOption(slv, name, v);
}
int AMPLSGetStrOption(AMPLS_MP_Solver* slv,
  const char* name,  const char** v) {
  std::string s;
  auto be = AMPLSGetBackend(slv);
  try {
    auto opt = be->GetOption(name);
    auto ii = ((AMPLS_MP__internal*)(slv->internal_info_));
    opt->GetValue(ii->optionvalue_);
    *v = ii->optionvalue_.data();
    return 0;
  }
  catch (const mp::OptionError&) {
    return 1;
  }
  catch (const std::exception&) {
    return 2;
  }
}
