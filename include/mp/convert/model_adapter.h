#ifndef MODEL_ADAPTER_H
#define MODEL_ADAPTER_H

#include "mp/convert/converter_query.h"

namespace mp {

/// Present an interface to "original" model
/// corresponding to the NL file
/// This is a workaround to use MP output facilities
template <class Model>
class ModelAdapter : public Model {
  int n_vars_orig_ = -1;
  int n_alg_cons_orig_= -1;
public:
  const Model& GetModel() const { return *this; }
  Model& GetModel() { return *this; }

  /// MP solution writing uses the following model size queries
  /** Returns the number of variables. */
  int num_vars() const { return n_vars_orig_; }
  /** Returns the number of algebraic constraints. */
  int num_algebraic_cons() const { return n_alg_cons_orig_; }
  void set_num_vars(int n) { n_vars_orig_ = n; }
  void set_num_alg_cons(int n) { n_alg_cons_orig_ = n; }


  ////////////////////////// SUFFIX I/O //////////////////////////////
  template <class T>
  void ReportSuffix(const SuffixDef<T>& sufdef,
                    ArrayRef<T> values) {
    if (values.empty())
      return;
    auto suf = FindOrCreateSuffix(sufdef);
    auto suf_size = suf.num_values();
    /// Check this because Converter or solver can add more variables
    assert(suf_size <= (int)values.size());
    for (auto i=suf_size; i--; ) {
      suf.set_value(i, values[i]);
    }
  }

  template <class T>
  ArrayRef<T> ReadSuffix(const SuffixDef<T>& sufdef) {
    auto suf = FindSuffix(sufdef);
    if (!suf)
      return {};
    return suf.get_values();
  }

  template <class T>
  BasicMutSuffix<T> FindSuffix(const SuffixDef<T>& sufdef) {
    auto main_kind = (suf::Kind)(sufdef.kind() & suf::KIND_MASK);
    auto suf_raw = Model::suffixes(main_kind).Find(sufdef.name());
    if (suf_raw)
      return Cast< BasicMutSuffix<T> >( suf_raw );
    return BasicMutSuffix<T>();
  }

  template <class T>
  BasicMutSuffix<T> FindOrCreateSuffix(const SuffixDef<T>& sufdef) {
    auto main_kind = (suf::Kind)(sufdef.kind() & suf::KIND_MASK);
    auto suf_raw = FindSuffix(sufdef);
    auto suf_size = GetSuffixSize(main_kind);    // can be < values.size()
    if (suf_raw) {
      suf_raw.or_kind(suf::OUTPUT);
      return suf_raw;
    }
    return Model::suffixes(main_kind).template
            Add<T>(sufdef.name(), sufdef.kind() | suf::OUTPUT,
                   suf_size, sufdef.table());
  }

  int GetSuffixSize(suf::Kind kind) {
    std::size_t size = 0;
    switch (kind) {
    default:
      MP_ASSERT(false, "invalid suffix kind");
      // Fall through.
    case suf::VAR:
      size = num_vars();
      break;
    case suf::CON:
      size = num_algebraic_cons();
      break;
    case suf::OBJ:
      size = Model::num_objs();
      break;
    case suf::PROBLEM:
      size = 1;
      break;
    }
    return static_cast<int>(size);
  }


};

} // namespace mp

#endif // MODEL_ADAPTER_H
