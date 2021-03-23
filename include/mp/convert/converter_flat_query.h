#ifndef CONVERTER_FLAT_QUERY_H
#define CONVERTER_FLAT_QUERY_H

/**
  * This is an implementation of a 'query class' for flat converters,
  * which provides solution handling, options, suffixes, etc.
  * To be used by a backend
  */

#include "mp/convert/converter_query.h"

namespace mp {

template <class Converter>
class FlatConverterQuery : public ConverterQuery {
  Converter& cvt_;
public:
  FlatConverterQuery(Converter& c) : cvt_(c) { }
  const Converter& GetCvt() const { return cvt_; }
  Converter& GetCvt() { return cvt_; }

  using Model = typename Converter::ModelType;

  const Model& GetPB() const { return GetCvt().GetModel(); }
  Model& GetPB() { return GetCvt().GetModel(); }

  IntSuffixHandler AddIntSuffix(fmt::StringRef name, int kind, int) override {
    return GetPB().AddIntSuffix(name, kind);
  }

};

} // namespace mp

#endif // CONVERTER_FLAT_QUERY_H
