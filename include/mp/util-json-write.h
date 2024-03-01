/*
 JSON write utils.

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

#ifndef UTILJSONWRITE_H
#define UTILJSONWRITE_H

#include <vector>
#include <string>
#include <tuple>
#include <type_traits>

namespace mp {

/// A lightweight JSON writer.
///
/// Does not require an intermediate representation of the data.
/// Similar to https://github.com/giacomodrago/minijson_writer
/// but the string buffer writer is a template parameter
/// with method write() conformant to fmt::format.
///
/// @tparam Formatter: consumes provided data (e.g., numbers/strings)
///   to produce a JSON string/file.
template <class Formatter>
class MiniJSONWriter {
public:
  /// Construct with provided Formatter \a wrt.
  MiniJSONWriter(Formatter& wrt) : wrt_(wrt) { }

  /// Declare Node.
  /// A node is a JSON tree node (value, array, dict.)
  using Node = MiniJSONWriter;

  /// Prefix++: make/ensure *this an array,
  /// add and return a new element.
  /// Useful when writing complex elements
  /// or non-supported types manually.
  Node operator++();

  /// operator[]: make/ensure *this a dictionary,
  /// add and return a new element at \a key.
  /// Allows the syntax `node[key] = val;`,
  /// or `node[key] << val1 << val2 ...;`,
  /// but also writing a complex subtree manually.
  Node operator[](std::string_view key);

  /// operator<<: make/ensure *this an array
  /// and write \a val a new element.
  ///
  /// Equivalent: `++(*this) = val; return *this;`
  ///
  /// @param val: scalar, tuple or container.
  ///   For non-supported types,
  ///   write elements manually into
  ///     `auto jw_child = ++jw;`,
  ///   or define, e.g.,
  ///     `Serialize(MiniJSONWriter , const YourType& )`,
  ///   to be used as follows:
  ///     `Serialize(++jw, obj[5]);`
  ///     (instead of `jw << obj[5]`.)
  ///
  /// @return *this.
  template <class Value>
  Node& operator<<(const Value& val)
  { ++(*this) = val; return *this; }

  /// operator=: write \a val as a whole.
  ///
  /// @param val: scalar, tuple or container.
  ///   For non-supported types,
  ///   write elements manually into
  ///     `auto jw_child = jw["data"];`,
  ///   or define, e.g.,
  ///     `Serialize(MiniJSONWriter , const YourType& )`,
  ///   to be used as follows:
  ///     `Serialize(jw["data"], obj_data);`
  ///     (instead of `jw["data"] = obj_data`.)
  ///
  /// @note Can be called only once on a single node.
  template <class Value>
  void operator=(const Value& val)
  { EnsureUnset(); Write(val); Close(); }

  /// Write a sequence between two iterators
  /// @return *this
  template <class It>
  Node& WriteSequence(It b, It e) {
    for ( ; b!=e; ++b)
      (*this) << (*b);
    return *this;
  }

  /// Close node.
  /// Call this before writing to any non-child nodes,
  /// if not exiting the scope which calls the destructor.
  void Close();

  /// Destructor.
  /// Closes node by RAII.
  ~MiniJSONWriter() { Close(); }

  /// Move construct
  MiniJSONWriter(MiniJSONWriter&& other)
    : wrt_(other.wrt_),
      kind_(other.kind_), n_written_(other.n_written_)
  { other.kind_ = Kind::Closed; }

protected:
  /// Construct a child.
  MiniJSONWriter(MiniJSONWriter& parent)
  : wrt_(parent.wrt_) { }

  /// Generic value write.
  template <class Value>
  void Write(const Value& val)
  { EnsureCanWrite(); DoWrite(val); }

  /// Make sure this node has not been written into.
  void EnsureUnset();
  /// If unset, mark scalar node.
  void MakeScalarIfUnset();
  /// Make sure this is an array node.
  void EnsureArray();
  /// Make sure this is a dictionary node.
  void EnsureDictionary();
  /// Make sure we can write another value/element.
  void EnsureCanWrite();

  /// Kind enum
  enum class Kind {
    Unset,
    Scalar,
    Array,
    Dict,
    Closed
  };

  /// Insert element separator if needed
  void InsertElementSeparator();

  void DoWrite(const char* s) { DoWriteString(s); }
  void DoWrite(const std::string& s) { DoWriteString(s); }

  template <typename Arithmetic,
            typename
            std::enable_if_t<std::is_arithmetic_v<Arithmetic>,
                             int> = 0> // C++17
  void DoWrite(Arithmetic v) { DoWriteScalar(v); }

  /// Write a container
  template <typename C,
      typename T = std::decay_t<
          decltype(*begin(std::declval<C>()))> >
  void DoWrite(const C& c) { WriteSequence(c.begin(), c.end()); }

  /// https://www.cppstories.com/2022/tuple-iteration-apply/
  template <typename... Arg>
  void DoWrite(const std::tuple<Arg...>& tup) {
    std::apply([this](const auto&... tupleArgs) {
      auto printElem = [this](const auto& x) {
        (*this) << x;
      };
      (printElem(tupleArgs), ...);
    }, tup
    );
  }

  template <class Value>
  void DoWriteScalar(const Value& val) {
    MakeScalarIfUnset();
    wrt_.write("{}", val);
    ++n_written_;
  }

  template <class Str>
  void DoWriteString(const Str& val) {
    MakeScalarIfUnset();
    wrt_.write("\"{}\"", val);
    ++n_written_;
  }

private:
  Formatter& wrt_;
  Kind kind_{Kind::Unset};
  int n_written_{};
};

}  // namespace mp

#endif // UTILJSONWRITE_H
