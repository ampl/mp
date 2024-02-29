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
  /// but also writing a complex subtree manually.
  Node operator[](const char* key);

  /// operator<<: make/ensure *this an array
  /// and write \a val a new element.
  ///
  /// @param val: scalar or container.
  ///   For non-supported types, define global method
  ///   Serialize(MiniJSONWriter&, const YourType& ).
  ///
  /// @return *this.
  template <class Value>
  Node& operator<<(const Value& val)
  { EnsureArray(); Write(val); return *this; }

  /// operator=: write \a val as a whole.
  ///
  /// @param val: scalar or container.
  ///   For non-supported types, define global method
  ///   Serialize(MiniJSONWrt& , const YourType& ).
  ///
  /// @note Can be called only once on a single node.
  template <class Value>
  void operator=(const Value& val)
  { EnsureUnset(); Write(val); Close(); }

  /// Write a sequence between two iterators
  template <class It>
  void WriteSequence(It b, It e) {
    for ( ; b!=e; ++b)
      (*this) << (*b);
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
  /// Make sure this is a scalar node.
  void EnsureScalar();
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

  void DoWrite(const char* s) { DoWriteScalar(s); }
  void DoWrite(const std::string& s) { DoWriteScalar(s); }

  template <typename Arithmetic,
            typename
            std::enable_if_t<std::is_arithmetic_v<Arithmetic>,
                             int> = 0> // C++17
  void DoWrite(Arithmetic v) { DoWriteScalar(v); }

  template <typename C,
      typename T = std::decay_t<
          decltype(*begin(std::declval<C>()))> >
  void DoWrite(const C& c) { WriteSequence(c.begin(), c.end()); }

  template <class Value>
  void DoWriteScalar(const Value& val)
  { EnsureScalar(); wrt_.write("{}", val); }

private:
  Formatter& wrt_;
  Kind kind_{Kind::Unset};
  int n_written_{};
};

}  // namespace mp

#endif // UTILJSONWRITE_H
