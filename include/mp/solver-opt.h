#ifndef SOLVEROPT_H
#define SOLVEROPT_H

#include <set>
#include <cstdint>
#include <cassert>

#include "mp/error.h"
#include "mp/format.h"


namespace mp {

/// Information about a possible option value.
struct OptionValueInfo {
  const char *value;
  const char *description;
  intptr_t data;  // Solver-specific data associated with this value.
};

/// A reference to an array of OptionValueInfo objects.
class ValueArrayRef {
 private:
  const OptionValueInfo *values_;
  int size_;

 public:
  ValueArrayRef() : values_(), size_() {}

  template <int SIZE>
  ValueArrayRef(const OptionValueInfo (&values)[SIZE], int offset = 0)
  : values_(values + offset), size_(SIZE - offset) {
    assert(offset >= 0 && offset < SIZE);
  }

  int size() const { return size_; }

  typedef const OptionValueInfo *iterator;

  iterator begin() const { return values_; }
  iterator end() const { return values_ + size_; }
};


namespace internal {

/// Formats the string s containing the reStructuredText (RST) markup
/// and writes it to w.
/// @param w:      stream writer used for output
/// @param indent: indentation to use for the formatted text
/// @param s:      string containing reStructuredText to format
/// @param values: information about possible option values to be formatted by the
///         value-table directive
void FormatRST(fmt::Writer &w, fmt::CStringRef s,
    int indent = 0, ValueArrayRef values = ValueArrayRef());

/// A helper class for implementing an option of type T.
template <typename T>
struct OptionHelper;

template <>
struct OptionHelper<int> {
  typedef int Arg;
  static void Write(fmt::Writer &w, Arg value) { w << value; }
  static int Parse(const char *&s);
  static int CastArg(fmt::LongLong value) { return static_cast<int>(value); }
};

template <>
struct OptionHelper<fmt::LongLong> {
  typedef fmt::LongLong Arg;
  static void Write(fmt::Writer &w, Arg value) { w << value; }
  static fmt::LongLong Parse(const char *&s) {
    return OptionHelper<int>::Parse(s); // TODO
  }
  static fmt::LongLong CastArg(fmt::LongLong value) { return value; }
};

template <>
struct OptionHelper<double> {
  typedef double Arg;
  static void Write(fmt::Writer &w, double value) { w << value; }
  static double Parse(const char *&s);
  static double CastArg(double value) { return value; }
};

template <>
struct OptionHelper<std::string> {
  typedef fmt::StringRef Arg;
  static void Write(fmt::Writer &w, const std::string &s) { w << s; }
  static std::string Parse(const char *&s);
  static fmt::StringRef CastArg(fmt::StringRef s) { return s; }
};

/// Raise OptionError
inline OptionError OptionTypeError(fmt::StringRef name, fmt::StringRef type) {
  return OptionError(
        fmt::format("Option \"{}\" is not of type \"{}\"", name, type));
}

}  // namespace internal


/// A solver option.
class SolverOption {
public:
  /// Constructs a SolverOption object.
  ///
  /// The solver option stores pointers to the passed name(s) and description and
  /// copies the strings for the names only. Normally the description is a
  /// string literal and has static storage duration but if this is not the
  /// case make sure that the string's lifetime is longer than that of the
  /// option object.
  ///
  /// The description should be written in a subset of reStructuredText (RST).
  /// Currently the following RST constructs are supported:
  ///
  /// * paragraphs
  /// * bullet lists
  /// * literal blocks
  /// * line blocks
  /// * the value-table directive (.. value-table::) which is replaced by a
  ///   table of option values as given by the values array
  ///
  /// names_list:  option names list
  ///   names can be wildcarded, e.g., "obj:*:method"
  ///   then setter/getter can use wildcard_ accessors
  /// description: option description
  /// values:      information about possible option values
  SolverOption(const char *names_list, const char *description,
      ValueArrayRef values = ValueArrayRef(), bool is_flag = false);

  virtual ~SolverOption() {}

  /// Return the option name.
  const char *name() const { return name_.c_str(); }

  /// Return the ASL (not qualified) name as the first
  /// inline synonym - or the name itself if no synonyms
  /// are defined
  const char* name_ASL() const;
  /// Returns the "inline" synonyms
  const std::vector<std::string>& inline_synonyms() const
  { return inline_synonyms_; }
  /// Add additional "inline" synonyms
  void add_synonyms_front(const char* names_list);
  void add_synonyms_back(const char* names_list);

  /// Wildcards
  bool is_wildcard() const { return wc_headtails_.size(); }
  /// Checks if matches, then saves key & body
  bool wc_match(const std::string& key);
  const std::string& wc_head() const
  { assert(is_wildcard()); return wc_headtails_[0].first; }
  const std::string& wc_tail() const
  { assert(is_wildcard()); return wc_headtails_[0].second; }
  const std::string& wc_key_last() const { return wc_key_last_; }
  const std::string& wc_keybody_last() const { return wc_body_last_; }
  /// Printing last parsed wc key in std form
  std::string wc_key_last__std_form() const
  { return wc_head() + wc_body_last_ + wc_tail(); }

  /// Return/set the option description.
  const char *description() const { return description_.c_str(); }
  void set_description(const char* d) { description_=d; }
  void add_to_description(const char* d) { description_ += d; }

  /// Returns the information about possible values.
  ValueArrayRef values() const { return values_; }

  /// Returns true if this option is a flag, i.e. it doesn't take a value.
  bool is_flag() const { return is_flag_; }

  /// Returns the option value.
  virtual void GetValue(fmt::LongLong &) const {
    throw internal::OptionTypeError(name_, "int");
  }
  virtual void GetValue(double &) const {
    throw internal::OptionTypeError(name_, "double");
  }
  virtual void GetValue(std::string &) const {
    throw internal::OptionTypeError(name_, "string");
  }

  virtual void GetValue(int &int_value) const {
    fmt::LongLong value = 0;
    GetValue(value);
    if (value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) {
      throw Error("Value {} doesn't fit in int", value);
    }
    int_value = static_cast<int>(value);
  }

  template <typename T>
  T GetValue() const {
    T value = T();
    GetValue(value);
    return value;
  }

  /// Sets the option value or throws InvalidOptionValue if the value is invalid.
  virtual void SetValue(fmt::LongLong) {
    throw internal::OptionTypeError(name_, "int");
  }
  virtual void SetValue(double) {
    throw internal::OptionTypeError(name_, "double");
  }
  virtual void SetValue(fmt::StringRef) {
    throw internal::OptionTypeError(name_, "string");
  }
  virtual void SetValue(int value) {
    fmt::LongLong long_value = value;
    SetValue(long_value);
  }

  /// Formats the option value. Throws OptionError in case of error.
  virtual void Write(fmt::Writer &w) = 0;

  /// Parses a string and sets the option value. Throws InvalidOptionValue
  /// if the value is invalid or OptionError in case of another error.
  virtual void Parse(const char *&s) = 0;

  virtual std::string echo() {
    if (is_wildcard())
      return wc_key_last__std_form();
    return name();
  }


private:
 std::string name_ {};
 std::vector<std::string> inline_synonyms_ {};
 std::string description_;

 /// Wildcard info
 /// Standard name's head/tail
 using WCHeadTail = std::pair<std::string, std::string>;
 std::vector<WCHeadTail> wc_headtails_;
 /// Last actual key parsed and it's '*' body
 std::string wc_key_last_, wc_body_last_;
 /// Assumes name constains '*'
 static WCHeadTail wc_split(const std::string& name);

 ValueArrayRef values_;
 bool is_flag_;
};


/// An exception thrown when an invalid value is provided for an option.
class InvalidOptionValue : public OptionError {
 private:
  template <typename T>
  static std::string Format(fmt::StringRef name, T value, fmt::StringRef msg) {
    if (0!=msg.size())
      return fmt::format("Invalid value \"{}\" for option \"{}\", {}",
                         value, name, msg);
    else
      return fmt::format("Invalid value \"{}\" for option \"{}\"", value, name);
  }

 public:
  template <typename T>
  InvalidOptionValue(fmt::StringRef name, T value, fmt::StringRef msg="")
  : OptionError(Format(name, value, msg)) {}

  template <typename T>
  InvalidOptionValue(const SolverOption &opt, T value, fmt::StringRef msg="")
  : OptionError(Format(opt.name(), value, msg)) {}
};

template <typename T>
class TypedSolverOption : public SolverOption {
 public:
  TypedSolverOption(const char *name, const char *description,
      ValueArrayRef values = ValueArrayRef())
  : SolverOption(name, description, values) {}

  void Write(fmt::Writer &w) { w << GetValue<T>(); }

  void Parse(const char *&s) {
    const char *start = s;
    T value = internal::OptionHelper<T>::Parse(s);
    if (*s && !std::isspace(*s)) {
      do ++s;
      while (*s && !std::isspace(*s));
      throw InvalidOptionValue(name(), std::string(start, s - start));
    }
    SetValue(value);
  }
};


/// A collection of solver options.
class SolverOptionManager {

public:

  /// Returns the number of options.
  int num_options() const { return static_cast<int>(options_.size()); }

  /// Returns the option with specified name.
  SolverOption *GetOption(const char *name) const {
    SolverOption *opt = FindOption(name);
    if (!opt)
      throw OptionError(fmt::format("Unknown option \"{}\"", name));
    return opt;
  }

  /// Handler should be a class derived from BasicSolver that will receive
  /// notifications about parsed options.
  template <typename Handler, typename T, typename AccessorT = T>
  class ConcreteOption : public TypedSolverOption<T> {
   private:
    typedef AccessorT (Handler::*Get)(const SolverOption &) const;
    typedef void (Handler::*Set)(
        const SolverOption &, typename internal::OptionHelper<AccessorT>::Arg);

    Handler &handler_;
    Get get_;
    Set set_;

   public:
    template <class Solver>
    ConcreteOption(const char *name, const char *description,
        Solver *s, Get get, Set set, ValueArrayRef values = ValueArrayRef())
    : TypedSolverOption<T>(name, description, values),
      handler_(static_cast<Handler&>(*s)), get_(get), set_(set) {}

    void GetValue(T &value) const { value = (handler_.*get_)(*this); }
    void SetValue(typename internal::OptionHelper<T>::Arg value) {
      (handler_.*set_)(*this,
                       internal::OptionHelper<AccessorT>::CastArg(value));
    }
  };
  template <typename Handler, typename T,
            typename Info, typename InfoArg = Info, typename AccessorT = T>
  class ConcreteOptionWithInfo : public TypedSolverOption<T> {
   private:
    typedef AccessorT (Handler::*Get)(const SolverOption &, InfoArg) const;
    typedef void (Handler::*Set)(
        const SolverOption &,
        typename internal::OptionHelper<AccessorT>::Arg, InfoArg);

    Handler &handler_;
    Get get_;
    Set set_;
    Info info_;

   public:
    template <class Solver>
    ConcreteOptionWithInfo(const char *name,
        const char *description, Solver *s, Get get, Set set, InfoArg info,
      ValueArrayRef values = ValueArrayRef())
    : TypedSolverOption<T>(name, description, values),
      handler_(static_cast<Handler&>(*s)), get_(get), set_(set), info_(info) {}
    ConcreteOptionWithInfo(const char *name, const char *description,
      Handler *s, Get get, Set set, InfoArg info, ValueArrayRef values = ValueArrayRef())
    : TypedSolverOption<T>(name, description, values),
      handler_(*s), get_(get), set_(set), info_(info) {}

    void GetValue(T &value) const { value = (handler_.*get_)(*this, info_); }
    void SetValue(typename internal::OptionHelper<T>::Arg value) {
      (handler_.*set_)(*this, internal::OptionHelper<AccessorT>::CastArg(value),
                       info_);
    }
  };

  /// Sets a text to be displayed before option descriptions.
  void set_option_header(const char *header) { option_header_ = header; }

  /// Add more text to be displayed before option descriptions.
  void add_to_option_header(const char *header_more) { option_header_ += header_more; }

  /// Smart ptr to an option
  typedef std::unique_ptr<SolverOption> OptionPtr;

  /// Add an option
  void AddOption(OptionPtr opt);

  /// Finds an option and returns a pointer to it if found or null otherwise.
  /// If wildcardvalues==true, wildcarded options should have values
  /// (i.e., parsing real input)
  SolverOption *FindOption(const char *name,
                           bool wildcardvalues=false) const;

  /// Returns the option header.
  const char *option_header() const { return option_header_.c_str(); }

  /// Flags for ParseOptions.
  enum {
    /// Don't echo options during parsing.
    NO_OPTION_ECHO = 1
  };

  /// Returns the value of an integer option.
  /// Throws OptionError if there is no such option or it has a different type.
  fmt::LongLong GetIntOption(const char *name) const {
    return GetOption(name)->GetValue<fmt::LongLong>();
  }

  // Sets the value of an integer option.
  // Throws OptionError if there is no such option or it has a different type.
  void SetIntOption(const char *name, fmt::LongLong value) {
    GetOption(name)->SetValue(value);
  }

  // Returns the value of a double option.
  // Throws OptionError if there is no such option or it has a different type.
  double GetDblOption(const char *name) const {
    return GetOption(name)->GetValue<double>();
  }

  // Sets the value of a double option.
  // Throws OptionError if there is no such option or it has a different type.
  void SetDblOption(const char *name, double value) {
    GetOption(name)->SetValue(value);
  }

  // Returns the value of a string option.
  // Throws OptionError if there is no such option or it has a different type.
  std::string GetStrOption(const char *name) const {
    return GetOption(name)->GetValue<std::string>();
  }

  // Sets the value of a string option.
  // Throws OptionError if there is no such option or it has a different type.
  void SetStrOption(const char *name, fmt::StringRef value) {
    GetOption(name)->SetValue(value);
  }


public:
  template <class Value>
  class StoredOption : public mp::TypedSolverOption<Value> {
    Value& value_;
  public:
    using value_type = Value;
    StoredOption(const char *name_list, const char *description,
                 Value& v, ValueArrayRef values = ValueArrayRef())
      : mp::TypedSolverOption<Value>(name_list, description, values), value_(v) {}

    void GetValue(Value &v) const override { v = value_; }
    void SetValue(typename internal::OptionHelper<Value>::Arg v) override
    { value_ = v; }
  };

  /// Simple stored option referencing a variable
  template <class Value>
  void AddStoredOption(const char *name, const char *description,
                       Value& value, ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(
                new StoredOption<Value>(
                  name, description, value, values)));
  }

  /// Simple stored option referencing a variable, min, max (TODO)
  template <class Value>
  void AddStoredOption(const char *name, const char *description,
                       Value& value, Value , Value ) {
    AddOption(OptionPtr(
                new StoredOption<Value>(
                  name, description, value, ValueArrayRef())));
  }

  /// Same: stored option referencing a variable, min, max (TODO)
  template <class Value>
  void AddOption(const char *name, const char *description,
                       Value& value, Value lb, Value ub) {
    AddStoredOption(name, description, value, lb, ub);
  }

  void ReplaceOptionDescription(const char* name, const char* desc) {
    auto pOption = FindOption(name);
    assert(pOption);
    pOption->set_description(desc);
  }

  void AddToOptionDescription(const char* name, const char* desc_add) {
    auto pOption = FindOption(name);
    assert(pOption);
    std::string to_add { "\n\n" };
    to_add += desc_add;
    pOption->add_to_description(to_add.c_str());
  }


  /// Add "inline" option synonyms
  /// The _Front version puts them in the front of the synonyms list
  /// and the 1st of them is used in the -a output for sorting
  void AddOptionSynonyms_Inline_Front(const char* names_list, const char* realName);
  void AddOptionSynonyms_Inline_Back(const char* names_list, const char* realName);

  /// Add an "out-of-line" synonym
  /// Creates extra entry under -=
  void AddOptionSynonyms_OutOfLine(const char* name, const char* realName);

  // Adds an integer option.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Int>
  void AddIntOption(const char *name,
    const char *description, Int (Handler::*get)(const SolverOption &) const,
      void (Handler::*set)(const SolverOption &, Int)) {
    AddOption(OptionPtr(new ConcreteOption<Handler, fmt::LongLong, Int>(
        name, description, this, get, set)));
  }

  // Adds an integer option with additional information.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Info>
  void AddIntOption(const char *name, const char *description,
      int (Handler::*get)(const SolverOption &, const Info &) const,
      void (Handler::*set)(const SolverOption &, int, const Info &),
      const Info &info) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<
                        Handler, fmt::LongLong, Info, const Info &, int>(
                          name, description, this, get, set, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Handler, typename Info>
  void AddIntOption(const char *name_list, const char *description,
      int (Handler::*get)(const SolverOption &, Info) const,
      void (Handler::*set)(const SolverOption &, int, Info), Info info) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<
                        Handler, fmt::LongLong, Info, Info, int>(
                          name_list, description, this, get, set, info)));
  }

  // Adds a double option.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler>
  void AddDblOption(const char *name, const char *description,
      double (Handler::*get)(const SolverOption &) const,
      void (Handler::*set)(const SolverOption &, double)) {
    AddOption(OptionPtr(new ConcreteOption<Handler, double>(
        name, description, this, get, set)));
  }

  // Adds a double option with additional information.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Info>
  void AddDblOption(const char *name, const char *description,
      double (Handler::*get)(const SolverOption &, const Info &) const,
      void (Handler::*set)(const SolverOption &, double, const Info &),
      const Info &info) {
    AddOption(OptionPtr(
        new ConcreteOptionWithInfo<Handler, double, Info, const Info &>(
            name, description, this, get, set, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Handler, typename Info>
  void AddDblOption(const char *name, const char *description,
      double (Handler::*get)(const SolverOption &, Info) const,
      void (Handler::*set)(const SolverOption &, double, Info), Info info) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<Handler, double, Info>(
            name, description, this, get, set, info)));
  }

  // Adds a string option.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler>
  void AddStrOption(const char *name, const char *description,
      std::string (Handler::*get)(const SolverOption &) const,
      void (Handler::*set)(const SolverOption &, fmt::StringRef),
      ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(new ConcreteOption<Handler, std::string>(
        name, description, this, get, set, values)));
  }

  // Adds a string option with additional information.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Info>
  void AddStrOption(const char *name, const char *description,
      std::string (Handler::*get)(const SolverOption &, const Info &) const,
      void (Handler::*set)(const SolverOption &, fmt::StringRef, const Info &),
      const Info &info, ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(
        new ConcreteOptionWithInfo<Handler, std::string, Info, const Info &>(
            name, description, this, get, set, info, values)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Handler, typename Info>
  void AddStrOption(const char *name, const char *description,
      std::string (Handler::*get)(const SolverOption &, Info) const,
      void (Handler::*set)(const SolverOption &, fmt::StringRef, Info),
      Info info, ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<Handler, std::string, Info>(
            name, description, this, get, set, info, values)));
  }

private:
  struct OptionNameLess {
    bool operator()(const SolverOption *lhs, const SolverOption *rhs) const;
  };

  std::string option_header_;
  typedef std::set<SolverOption*, OptionNameLess> OptionSet;
  OptionSet options_;


public:
  /// Option iterator.
  class option_iterator :
    public std::iterator<std::forward_iterator_tag, SolverOption> {
   private:
    OptionSet::const_iterator it_;

    friend class SolverOptionManager;

    explicit option_iterator(OptionSet::const_iterator it) : it_(it) {}

   public:
    option_iterator() {}

    const SolverOption &operator*() const { return **it_; }
    const SolverOption *operator->() const { return *it_; }

    option_iterator &operator++() {
      ++it_;
      return *this;
    }

    option_iterator operator++(int ) {
      option_iterator it(*this);
      ++it_;
      return it;
    }

    bool operator==(option_iterator other) const { return it_ == other.it_; }
    bool operator!=(option_iterator other) const { return it_ != other.it_; }
  };

  option_iterator option_begin() const {
    return option_iterator(options_.begin());
  }
  option_iterator option_end() const {
    return option_iterator(options_.end());
  }

  ~SolverOptionManager();
};

}  // namespace mp

#endif // SOLVEROPT_H
