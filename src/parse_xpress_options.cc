#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cfloat>
#include <cctype>

#include "wrap_xml.h"
#include "mp/utils-string.h"

#define RAISE(streammsg) do { \
  std::cerr << "ERROR: " << streammsg << std::endl; \
  std::exit(1); \
} while (0)

#define SET_ONCE(param, var, val) { \
  if ((var).size()) \
    RAISE("Param " << param << ": repeated data for " << #var); \
  if (std::string(val).empty()) \
    RAISE("Param " << param << ": empty value for " << #var); \
  (var) = (val); \
}

namespace mp {

/// Parameter
class Param {
public:
  std::string name_main_, name2_, type_;
  std::string xml_name_;      // name as in the XML, without key prefix
  std::string descr_;
  std::string topic_;
  std::string category_;
  std::string default_;
  std::string value_type_;    // num/bit
  std::map<std::string, std::string> values_map_;
  std::list< std::pair<std::string, std::string> > values_;
  /// check if ok
  std::string Check() const {
    if (name_main_.empty())
      return "No main name";
    if (descr_.empty())
      return "No description";
    if (topic_.empty())
      return "No topic";
    if (category_.empty())
      return "No category";
    if (default_.empty() && IsControl())
      return "No default value for a control";
    return "";
  }
  /// 1st name
  const char* NameMain() const { return name_main_.c_str(); }
  /// Is a control?
  bool IsControl() const { return "Control"==category_; }
};

/// Control types read from the Xpress C headers (xprs.h, xslp.h).
/// The XML parameter database is the source of names and documentation,
/// but its type attribute is not always right (e.g., OBJSCALEFACTOR is
/// declared as double in the 9.8 database but is an integer control).
/// The headers group the controls in sections such as
/// "/* Double control parameters */", which we use as the authority.
class HeaderControlTypes {
  std::map<std::string, std::string> types_;     // name -> int/double/string/attr
  std::map<std::string, std::string> aliases_;   // name -> other macro
  std::map<std::string, bool> deprecated_;
public:
  /// Read a header file. Returns false if it cannot be opened.
  bool ReadFile(const char* fln) {
    std::ifstream ifs(fln);
    if (!ifs)
      return false;
    std::string line, sec;           // sec: "", int, double, string, attr
    while (std::getline(ifs, line)) {
      auto cmt = line.find("/*");
      auto def = line.find("#define");
      if (std::string::npos != cmt && (std::string::npos==def || cmt<def)) {
        auto txt = MakeLower(line.substr(cmt));
        if (std::string::npos != txt.find("attribute")
            && std::string::npos == txt.find("control"))
          sec = "attr";
        else if (std::string::npos != txt.find("control")) {
          if (std::string::npos != txt.find("double"))
            sec = "double";
          else if (std::string::npos != txt.find("string"))
            sec = "string";
          else
            sec = "int";
        } else
          sec = "";                      // e.g. "/* Knitro's parameters */"
      }
      if (std::string::npos != def) {
        std::string rest = line.substr(def + 7);
        auto p1 = rest.find_first_not_of(" \t");
        if (std::string::npos == p1)
          continue;
        auto p2 = rest.find_first_of(" \t", p1);
        auto name = rest.substr(p1, p2-p1);
        if (0 != name.find("XPRS_") && 0 != name.find("XSLP_")
            && 0 != name.find("XKTR_"))
          continue;
        std::string val;
        if (std::string::npos != p2) {
          auto p3 = rest.find_first_not_of(" \t", p2);
          if (std::string::npos != p3) {
            auto p4 = rest.find_first_of(" \t", p3);
            val = rest.substr(p3, p4-p3);
          }
        }
        deprecated_[name] = (std::string::npos != line.find("Deprecated"));
        if (val.size() && (std::isdigit((unsigned char)val[0]) || '('==val[0]))
          types_[name] = sec;
        else if (val.size() && (0==val.find("XPRS_") || 0==val.find("XSLP_")))
          aliases_[name] = val;
      }
    }
    return true;
  }
  /// Header type of a control: int/double/string, "attr" for attributes,
  /// "" if unknown.
  std::string Type(const std::string& name) const {
    auto nm = name;
    for (int i=0; i<5; ++i) {              // follow aliases
      auto it = types_.find(nm);
      if (types_.end() != it)
        return it->second;
      auto ia = aliases_.find(nm);
      if (aliases_.end() == ia)
        break;
      nm = ia->second;
    }
    return "";
  }
  /// Is deprecated according to the header?
  bool IsDeprecated(const std::string& name) const {
    auto it = deprecated_.find(name);
    return deprecated_.end() != it && it->second;
  }
  /// Any headers read?
  bool Empty() const { return types_.empty(); }
  /// Lowercase
  static std::string MakeLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
  }
};


/// Param list translator
class ParamListTranslator {
  std::string key_prefix_;
  bool nonlinear_;                 // nonlinearparams.xml: slp/mislp/xktr prefixes
  const HeaderControlTypes& hdr_types_;

  std::ofstream ofs_;
  std::string filename_, xmlfile_, classname_, hdr_, ftr_;

  std::map<std::string, Param> params_;
  std::map<std::string, int> prefixes_;
public:
  /// Construct
  ParamListTranslator(const char* prefix, bool nonlinear,
                      const HeaderControlTypes& ht)
      : key_prefix_(prefix ? prefix : ""), nonlinear_(nonlinear),
        hdr_types_(ht) { }
  /// Start output
  bool Start(const char* xmlfile, const char* filename, const char* classname) {
    classname_ = classname;
    filename_ = filename;
    xmlfile_ = xmlfile;
    auto pos = xmlfile_.find_last_of("/\\");
    if (std::string::npos != pos)
      xmlfile_ = xmlfile_.substr(pos+1);
    if (OpenFile(filename)) {
      return true;
    }
    return false;
  }

  /// Destroy
  ~ParamListTranslator() {
    if (ofs_) {
      ofs_ << GetHeader(xmlfile_, classname_);
      WriteControlParams(ofs_);
      ofs_ << GetFooter(classname_);
      ofs_.close();
      OutputPrefixes();
    }
  }

  /// Add new parameter
  void AddParam(Param prm) {
    if (kControlsNotExposed.count(prm.name_main_)) {
      std::cerr << "NOTE: Skipping param " << prm.name_main_
                << ": not exposed as an option" << std::endl;
      return;
    }
    prm.xml_name_ = prm.name_main_;
    auto name_w_pref = key_prefix_;
    if (name_w_pref.size())
      name_w_pref += '_';
    name_w_pref += prm.NameMain();
    if (params_.end() != params_.find(name_w_pref)) {
      RAISE( "parameter '"
            << name_w_pref << "' repeated.");
    }
    prm.name_main_ = name_w_pref;
    auto ht = hdr_types_.Type(name_w_pref);
    if ("attr" == ht) {
      std::cerr << "WARNING: Skipping param " << name_w_pref
                << ": the header declares it as an attribute" << std::endl;
      return;
    }
    if (hdr_types_.IsDeprecated(name_w_pref)) {
      std::cerr << "WARNING: Skipping param " << name_w_pref
                << ": deprecated in the header" << std::endl;
      return;
    }
    if (ht.size()) {
      std::string xt = ("integer"==prm.type_) ? "int" : prm.type_;
      if (xt != ht) {
        std::cerr << "WARNING: Param " << name_w_pref
                  << ": XML type '" << prm.type_
                  << "' but the header section says '" << ht
                  << "'; using the header type" << std::endl;
        prm.type_ = ("int"==ht) ? "integer" : ht;
      }
    }
    params_[name_w_pref] = std::move(prm);
  }

protected:
  /// Open output file
  bool OpenFile(const char* name) {
    ofs_.open(name);
    if (!ofs_) {
      std::cerr << "Opening '" << name
                << "' for output failed: " << std::endl;
      std::perror("");
      return false;
    }
    return true;
  }
  /// Output file header part
  const char* GetHeader(
      const std::string& filename,
      const std::string& classname) {
    hdr_ =
        "#include <climits>\n"
        "#include <cfloat>\n"
          "\n"
        "#include \"mp/common.h\"\n"
        "#include \"mp/error.h\"\n"
        "#include \"mp/backend-std.h\"\n"
           "\n"
        "extern \"C\" {\n"
        "  #include \"xprs.h\"\n"
        "  #include \"xslp.h\"\n"
        "}\n"
           "\n\n"
        "namespace mp {\n"
        "\n"
        "/// A mix-in class to add Xpress parameters.\n"
        "/// Generated by parse-xpress-options from '" + filename +
        "'\n/// (parameter types verified against the Xpress headers).\n"
        "/// Do not edit; regenerate when the Xpress SDK changes.\n"
        "///\n"
        "template <class Impl>\n"
        "class Compiled" + classname + "Options {\n";
    return hdr_.c_str();
  }
  /// Output file footer part
  const char* GetFooter(
      const std::string& classname) {
    ftr_ =
        "  }  // Add" + classname + "Options()\n\n"
                      "};  // class Compiled" + classname +
        "Options\n\n"
        "}  // namespace mp\n";
    return ftr_.c_str();
  }
  /// Write the control params
  void WriteControlParams(std::ostream& os) {
    os <<
        "public:\n"
        "  /// Add up to " << params_.size()
       << " '" << classname_
       << "' parameters\n"
          "  void Add" << classname_
       << "Options() {"
          "\n\n";

    // AddSolverOption() should handle duplicate solver-side params
    //   @todo and this for all solvers
    for (const auto& prmval: params_) {
      const auto& prm = prmval.second;
      os << "#ifdef " << prm.NameMain() << "\n";

      auto prefix = GetOptionPrefix(prm.topic_);
      ++prefixes_[prefix];
      os << "    MPD( AddSolverOption_MergeDuplicates(\""
         << MakeOptionNames(prm, prefix) << "\",\n"
            "      \"" << prm.descr_ << "\"\n";
      if (prm.values_.size()) {
        os << "      \"\\n\\n\"\n      \"Values (default: "
              << prm.default_ << "):\\n\"";
        for (const auto& valdescr: prm.values_) {
          if (valdescr.first.empty())          // group heading
            os << "\n      \"\\n\\n" << valdescr.second << "\\n\"";
          else
            os << "\n      \"\\n- (" << valdescr.first
               << ")  " << valdescr.second << "\"";
        }
      } else
        os << "      \"\\n\\nDefault: " << prm.default_ << "\"";
      os << ",\n"
         << "      ";
      os << prm.NameMain();
      if ("double" == prm.type_)
        os << ", -DBL_MAX, DBL_MAX";
      else if ("integer" == prm.type_)
        os << ", INT_MIN, INT_MAX";
      else if ("string" == prm.type_)
      {}
      else
        RAISE("unknown param type " << prm.type_
                                           << " for param "
                                           << prm.NameMain());
      os << ") );\n";

      os << "#endif  // ifdef " << prm.NameMain() << "\n\n";
    }
  }

  /// Option names for a control:
  /// "<prefix>:<short> <plain> <MACRO> [<name2>]", where <short> is the
  /// lowercase XML name without the XSLP_/XKTR_PARAM_ fragments and
  /// <plain> is the lowercase control name (the ASL keyword spelling
  /// for optimizer controls, the lowercase macro for nonlinear ones).
  std::string MakeOptionNames(const Param& prm, const std::string& prefix) const {
    auto lower = MakeOptionName(prm.xml_name_);
    auto plain = key_prefix_.size() ? lower : MakeOptionName(prm.NameMain());
    auto shrt = lower;
    for (const char* frag: {"xslp_", "xktr_param_"})
      if (0 == shrt.find(frag))
        shrt = shrt.substr(std::strlen(frag));
    auto alias = kShortNameOverrides.find(prm.xml_name_);
    if (kShortNameOverrides.end() != alias)
      shrt = alias->second;
    std::string names = prefix + ':' + shrt;
    if (shrt != plain)
      names += ' ' + shrt;
    names += ' ' + plain;
    auto legacy = kLegacyKeywords.find(prm.xml_name_);
    if (kLegacyKeywords.end() != legacy)
      names += ' ' + legacy->second;
    names += ' ' + prm.name_main_;
    if (prm.name2_.size())
      names += ' ' + prm.name2_;
    return names;
  }

  /// Controls present in the database that make no sense as options
  const std::set<std::string> kControlsNotExposed { "VERSION" };

  /// Short names differing from the XML name
  /// (established AMPL keywords of the ASL driver)
  const std::map<std::string, std::string> kShortNameOverrides {
    {"IFCHECKCONVEXITY", "convexitychk"},
  };

  /// Keywords of the ASL driver that differ from the control name,
  /// added as synonyms
  const std::map<std::string, std::string> kLegacyKeywords {
    {"EIGENVALUETOL", "eigenvaltol"},
    {"FEASIBILITYPUMP", "feaspump"},
    {"FEASTOLPERTURB", "feastol_perturb"},
    {"FEASTOLTARGET", "feastol_target"},
    {"HEURSEARCHEFFORT", "heureffort"},
    {"HEURSEARCHFREQ", "heursearch"},
    {"HEURSEARCHROOTSELECT", "heurroot"},
    {"HEURSEARCHTREESELECT", "heurtree"},
    {"HISTORYCOSTS", "pseudocost_ud"},
    {"LPREFINEITERLIMIT", "lpref_itlim"},
    {"MAXMEMORYSOFT", "maxmemory"},
    {"MAXSCALEFACTOR", "maxlogscale"},
    {"MIPCONCURRENTNODES", "mipconcurnodes"},
    {"MIPCONCURRENTSOLVES", "mipconcursolves"},
    {"MIPRESTARTGAPTHRESHOLD", "miprestartgaptol"},
    {"MIPTERMINATIONMETHOD", "mipstop"},
    {"OPTIMALITYTOLTARGET", "opttol_target"},
    {"PRECONEDECOMP", "conedecomp"},
    {"PRECONVERTSEPARABLE", "preconvertsep"},
    {"QSIMPLEXOPS", "mipops"},
    {"TREECOMPRESSION", "treecompress"},
    {"TREECUTSELECT", "treecuts"},
    {"TREEDIAGNOSTICS", "treeoutlev"},
    {"TREEMEMORYLIMIT", "treememlimit"},
    {"TREEMEMORYSAVINGTARGET", "treememtarget"},
  };

  std::string MakeOptionName(std::string s) const {
    std::transform(s.begin(), s.end(), s.begin(),
                   // static_cast<int(*)(int)>(std::tolower)         // wrong
                   // [](int c){ return std::tolower(c); }           // wrong
                   // [](char c){ return std::tolower(c); }          // wrong
                   [](unsigned char c){ return std::tolower(c); } // correct
                   );
    return s;
  }

  /// MP option prefix for an XML topic string such as "LP, Barrier".
  /// The first matching entry of the table wins, so specific topics
  /// (Global, PDHG, Quadratic, Barrier, Cuts, ...) precede generic ones
  /// (Presolve, LP, Branch and Bound). Unmatched topics (Misc, Logging,
  /// System, Memory, Parallel, Callback, File IO, ...) go to "tech".
  /// Nonlinear controls use the engine prefixes xktr, mislp and slp only.
  std::string GetOptionPrefix(const std::string& topic) const {
    using Table = std::vector< std::pair<const char*, const char*> >;
    static const Table kNonlinear {
      {"Knitro", "xktr"},
      {"MISLP", "mislp"},
    };
    static const Table kOptimizer {
      {"Global", "global"},
      {"Primal Dual Hybrid Gradient", "pdhg"},
      {"Multiobjective", "obj:multi"},
      {"Infeasibility", "iis"},
      {"Tuner", "tech"},
      {"Quadratic", "qp"},
      {"Cuts", "cut"},
      {"Barrier", "bar"},
      {"Heuristics", "mip"},
      {"Branching", "mip"},
      {"Branch and Bound", "mip"},
      {"MIP Entities", "mip"},
      {"Root", "mip"},
      {"Presolve", "pre"},
      {"Propagation", "pre"},
      {"Problem Transformation", "pre"},
      {"Piecewise Linear", "pre"},
      {"Simplex", "lp"},
      {"LP", "lp"},
      {"Limits", "lim"},
      {"Solution Process", "alg"},
      {"Solution Refinement", "alg"},
      {"Numerics", "alg"},
      {"Tolerances", "alg"},
      {"Solution", "sol"},
    };
    const Table& table = nonlinear_ ? kNonlinear : kOptimizer;
    for (const auto& entry: table)
      if (std::string::npos != topic.find(entry.first))
        return entry.second;
    return nonlinear_ ? "slp" : "tech";
  }

  /// Print prefixes used
  void OutputPrefixes() const {
    std::cout << "Prefixes used:\n";
    for (const auto& p: prefixes_)
      std::cout << "  " << p.first << ":\t" << p.second << std::endl;
  }
};


/// Remove HTML tags such as <b>...</b>
std::string StripTags(const std::string& s) {
  std::string r;
  bool in_tag = false;
  for (char c: s) {
    if ('<' == c) in_tag = true;
    else if ('>' == c) in_tag = false;
    else if (!in_tag) r += c;
  }
  return r;
}

/// Collapse runs of blanks, trim both ends
std::string CollapseBlanks(const std::string& s) {
  std::string r;
  for (char c: s)
    if (' ' != c || (r.size() && ' ' != r.back()))
      r += c;
  while (r.size() && ' ' == r.back())
    r.pop_back();
  return r;
}

/// Replace line breaks
std::string GetText(BasicTreeWalker& wlk) {
  auto txt = wlk.GetText();
  // Replace all '\n' characters with spaces
  std::replace(txt.begin(), txt.end(), '\n', ' ');
  // Replace all '\r' characters with spaces
  std::replace(txt.begin(), txt.end(), '\r', ' ');
  std::replace(txt.begin(), txt.end(), '"', '\'');   //  " -> '
  txt = StripTags(txt);
  return CollapseBlanks(txt);
}


/// Specialize for value list
class ValueListHandler
    : public BasicNodeHandler {
  Param& prm_;
public:
  /// Construct
  ValueListHandler(Param& prm) : prm_(prm) { }

  /// Handle atribute
  void HandleAttribute(
      const char* name, const char* value) override {
    if (0 == std::strcmp("type", name)) {
      prm_.value_type_ = value;
    } else {
      RAISE("unknown <paramValues>'s attr: " << name);
    }
  }

  /// Handle param's entry
  void HandleSubnode(
      const char* name, BasicTreeWalker& wlk) override {
    if (0 == std::strcmp("paramVal", name)) {
      auto vstr = wlk.GetAttribute("value");
      if (!vstr)
        RAISE("For param " << prm_.NameMain()
                           << ", a value list entry has no value attr");
      if ('\0' == *vstr) {                   // a group heading
        prm_.values_.push_back({"", GetText(wlk)});
        return;
      }
      if (prm_.values_map_.end() != prm_.values_map_.find(vstr))
        std::cerr << "WARNING: For param " << prm_.NameMain()
                           << ", value list entry '"
                  << vstr << "' is repeated" << std::endl;
      auto vtext = GetText(wlk);
      if (vtext.empty()) {
        std::cerr << "WARNING: For param " << prm_.NameMain()
          << ", value list entry '"
          << vstr << "' has no description"
          << std::endl;
      }
      prm_.values_map_[vstr] = vtext;
      prm_.values_.push_back({vstr, vtext});
    } else {
      RAISE("unknown <paramValues>'s subnode: " << name);
    }
  }
};


/// Specialize for entries of <param>
class ParamHandler
    : public BasicNodeHandler {
  Param prm_;
public:
  /// Handle atribute
  void HandleAttribute(
      const char* name, const char* value) override {
    if (0 == std::strcmp("name", name)) {
      SET_ONCE(name, prm_.name_main_, value);
    } else if (0 == std::strcmp("name2", name)) {
      SET_ONCE(prm_.NameMain(), prm_.name2_, value);
    } else if (0 == std::strcmp("type", name)) {
      SET_ONCE(prm_.NameMain(), prm_.type_, value);
    } else  {
      RAISE("unknown <param>'s attribute: " << name);
    }
  }

  /// Handle param's entry
  void HandleSubnode(
      const char* name, BasicTreeWalker& wlk) override {
    if (0 == std::strcmp("paramDescr", name)) {
      SET_ONCE(prm_.NameMain(), prm_.descr_, GetText(wlk));
    } else if (0 == std::strcmp("paramTopic", name)) {
      SET_ONCE(prm_.NameMain(), prm_.topic_, GetText(wlk));
    } else if (0 == std::strcmp("paramCategory", name)) {
      SET_ONCE(prm_.NameMain(), prm_.category_, GetText(wlk));
    } else if (0 == std::strcmp("paramDefault", name)) {
      SET_ONCE(prm_.NameMain(), prm_.default_, GetText(wlk));
    } else if (0 == std::strcmp("paramValues", name)) {
      ValueListHandler vlh(prm_);
      wlk.Walk(vlh);
    } else if (0 == std::strcmp("paramNote", name)) {
      // skip
    } else {
      RAISE("unknown <param>'s subnode: " << name);
    }
  }

  /// const Param&
  const Param& ParamRef() const { return prm_; }

  /// Move out the param
  Param&& MoveOutParam() { return std::move(prm_); }
};


/// Specialize for <param>s, entries of <paramList>
class ParamListHandler
    : public BasicNodeHandler {
  ParamListTranslator& plt_;
public:
  /// Construct
  ParamListHandler(ParamListTranslator& plt) : plt_(plt) { }
  /// Handle param
  void HandleSubnode(
      const char* name, BasicTreeWalker& wlk) override {
    if (0 == std::strcmp("param", name)) {
      ParamHandler prmh;
      wlk.Walk(prmh);
      if (prmh.ParamRef().IsControl()) {     // Control param only
        auto err = prmh.ParamRef().Check();
        if (err.size())
          std::cerr << "WARNING: Skipping param "
                    << prmh.ParamRef().NameMain()
                    << ": " << err << std::endl;
        else
          plt_.AddParam(prmh.MoveOutParam());
      }
    } else {
      RAISE("unknown <paramList>'s subnode: " << name);
    }
  }
};

}  // namespace mp


int main(int argc, const char** argv) {
  if (argc<4)
    RAISE( "Usage: parse-xpress-options <params.xml> <output.h> <ClassName>"
          " [<key prefix>|-] [--nonlinear] [<xprs.h> [<xslp.h> ...]]\n"
          "  <params.xml>  XML parameter database from the Xpress SDK\n"
          "  <output.h>    C++ header to write (mix-in class Compiled<ClassName>Options)\n"
          "  <key prefix>  common macro prefix, e.g. XPRS for optimizerparams.xml;"
          " '-' for none\n"
          "  --nonlinear   nonlinearparams.xml: use the slp/mislp/xktr prefixes\n"
          "  <xprs.h> ...  Xpress C headers; when given, control types are taken"
          " from them\n"
          "                and controls they declare as attributes or deprecated"
          " are skipped." );

  std::cout << "Processing XML parameter database file '"
            << argv[1] << "' ..." << std::endl;

  const char* prefix = (argc>4 && std::strcmp("-", argv[4])) ? argv[4] : nullptr;
  bool nonlinear = false;
  mp::HeaderControlTypes hdr_types;
  for (int i=5; i<argc; ++i) {
    if (0 == std::strcmp("--nonlinear", argv[i])) {
      nonlinear = true;
    } else if (!hdr_types.ReadFile(argv[i]))
      RAISE("Cannot read header file '" << argv[i] << "'");
    else
      std::cout << "   Control types read from '" << argv[i] << "'" << std::endl;
  }

  auto pwlk = mp::MakeDefaultXMLWalker();
  if (pwlk->ReadFile(argv[1])) {
    if ( std::strcmp("paramList", pwlk->GetName()) )
      RAISE("Unknown top-level entry: " << pwlk->GetName());

    mp::ParamListTranslator plt(prefix, nonlinear, hdr_types);
    if (plt.Start(argv[1], argv[2], argv[3])) {

      mp::ParamListHandler prmlh(plt);  // to handle individual <param>s
      pwlk->Walk(prmlh);

      std::cout << "Writing header file '"
                << argv[2] << "' with mix-in class '"
                << argv[3] << "' ..." << std::endl;
      if (prefix)
        std::cout << "   (Key name prefix '" << prefix
                  << "' applied)" << std::endl;
    }
  } else
    RAISE("Error reading input file");

  std::cout << "Done." << std::endl;
  return 0;
}
