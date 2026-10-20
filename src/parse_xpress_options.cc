#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>
#include <map>
#include <list>
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
  const HeaderControlTypes& hdr_types_;

  std::ofstream ofs_;
  std::string filename_, xmlfile_, classname_, hdr_, ftr_;

  std::map<std::string, Param> params_;
  std::map<std::string, int> prefixes_;
public:
  /// Construct
  ParamListTranslator(const char* prefix, const HeaderControlTypes& ht)
      : key_prefix_(prefix ? prefix : ""), hdr_types_(ht) { }
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
      os << "    MPD( AddSolverOption_MergeDuplicates(\"";
      auto nm1 = MakeOptionName(prm.NameMain());
      if (0 == nm1.find(prefix)) {
        nm1 = nm1.substr(prefix.size()         // bar:alg
                         + (nm1.size()>prefix.size()
                            && '_'==nm1[prefix.size()])); // xktr_param..
      }
      os << prefix << ':' << nm1;
      os << ' ' << (prm.NameMain());           // no prefix
      if (prm.name2_.size())
        os << ' ' << (prm.name2_);             // no prefix
      os << "\",\n"
            "      \"" << prm.descr_ << "\"\n";
      if (prm.values_.size()) {
        os << "      \"\\n\\n\"\n      \"Values (default: "
              << prm.default_ << "):\\n\"";
        for (const auto& valdescr: prm.values_) {
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

  std::string MakeOptionName(std::string s) const {
    std::transform(s.begin(), s.end(), s.begin(),
                   // static_cast<int(*)(int)>(std::tolower)         // wrong
                   // [](int c){ return std::tolower(c); }           // wrong
                   // [](char c){ return std::tolower(c); }          // wrong
                   [](unsigned char c){ return std::tolower(c); } // correct
                   );
    return s;
  }

  std::string GetOptionPrefix(const std::string& topic) const {
    if (std::string::npos != topic.find("Global")) {
      return "global";
    }
    if (std::string::npos != topic.find("nitro")) {
      return "xktr";
    }
    if (std::string::npos != topic.find("Cuts")) {
      return "cut";
    }
    if (std::string::npos != topic.find("euristic")) {
      return "heur";
    }
    if (std::string::npos != topic.find("arrier")) {
      return "bar";
    }
    if (std::string::npos != topic.find("resolve")) {
      return "pre";
    }
    if (std::string::npos != topic.find("ropagation")) {
      return "pre";
    }
    if (std::string::npos != topic.find("unction")) {
      return "func";
    }
    if (std::string::npos != topic.find("uadrat")) {
      return "qp";
    }
    if (std::string::npos != topic.find("erivat")) {
      return "diff";
    }
    if (begins_with(topic, "MISLP")) {
      return "mislp";
    }
    if (std::string::npos != topic.find("Misc")) {
      return "tech";
    }
    if (std::string::npos != topic.find("ranch")) {
      return "mip";
    }
    if (std::string::npos != topic.find("arallel")) {
      return "tech";
    }
    if (std::string::npos != topic.find("eterminism")) {
      return "tech";
    }
    if (std::string::npos != topic.find("roblem Creation")) {
      return "prob";
    }
    if (std::string::npos != topic.find("File IO")) {
      return "tech";
    }
    if (std::string::npos != topic.find("rocess")) {
      return "alg";
    }
    if (std::string::npos != topic.find("ultistart")) {
      return "alg";
    }
    if (std::string::npos != topic.find("ultiobj")) {
      return "obj:multi";
    }
    if (std::string::npos != topic.find("rimal Dual")) {
      return "pdhg";
    }
    if (std::string::npos != topic.find("allback")) {
      return "tech";
    }
    if (std::string::npos != topic.find("ompute")) {
      return "tech";
    }
    if (std::string::npos != topic.find("Tuner")) {
      return "tech";
    }
    if (std::string::npos != topic.find("emory")) {
      return "tech";
    }
    if (std::string::npos != topic.find("SLP")) {
      return "slp";
    }
    if (std::string::npos != topic.find("LP")) {
      return "lp";
    }
    if (std::string::npos != topic.find("Root")) {
      return "mip";
    }
    auto p1 = topic.find_first_of(", ");
    auto p_end = (std::string::npos!=p1) ? p1 : topic.size();
    if (p_end>3)
      p_end=3;
    return MakeOptionName(topic.substr(0, p_end));
  }

  /// Print prefixes used
  void OutputPrefixes() const {
    std::cout << "Prefixes used:\n";
    for (const auto& p: prefixes_)
      std::cout << "  " << p.first << ":\t" << p.second << std::endl;
  }
};


/// Replace line breaks
std::string GetText(BasicTreeWalker& wlk) {
  auto txt = wlk.GetText();
  // Replace all '\n' characters with spaces
  std::replace(txt.begin(), txt.end(), '\n', ' ');
  // Replace all '\r' characters with spaces
  std::replace(txt.begin(), txt.end(), '\r', ' ');
  std::replace(txt.begin(), txt.end(), '"', '\'');   //  " -> '
  return txt;
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
      if ('\0' == *vstr) {
        std::cerr << "WARNING: For param " << prm_.NameMain()
                  << ", a value list entry has empty value" << std::endl;
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
          " [<key prefix>|-] [<xprs.h> [<xslp.h> ...]]\n"
          "  <params.xml>  XML parameter database from the Xpress SDK\n"
          "  <output.h>    C++ header to write (mix-in class Compiled<ClassName>Options)\n"
          "  <key prefix>  common macro prefix, e.g. XPRS for optimizerparams.xml;"
          " '-' for none\n"
          "  <xprs.h> ...  Xpress C headers; when given, control types are taken"
          " from them\n"
          "                and controls they declare as attributes or deprecated"
          " are skipped." );

  std::cout << "Processing XML parameter database file '"
            << argv[1] << "' ..." << std::endl;

  const char* prefix = (argc>4 && std::strcmp("-", argv[4])) ? argv[4] : nullptr;
  mp::HeaderControlTypes hdr_types;
  for (int i=5; i<argc; ++i) {
    if (!hdr_types.ReadFile(argv[i]))
      RAISE("Cannot read header file '" << argv[i] << "'");
    std::cout << "   Control types read from '" << argv[i] << "'" << std::endl;
  }

  auto pwlk = mp::MakeDefaultXMLWalker();
  if (pwlk->ReadFile(argv[1])) {
    if ( std::strcmp("paramList", pwlk->GetName()) )
      RAISE("Unknown top-level entry: " << pwlk->GetName());

    mp::ParamListTranslator plt(prefix, hdr_types);
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
