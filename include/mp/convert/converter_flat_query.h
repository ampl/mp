#ifndef CONVERTER_FLAT_QUERY_H
#define CONVERTER_FLAT_QUERY_H

/**
  * This is an implementation of a 'query class' for flat converters,
  * which provides solution handling, suffixes, etc.
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

  using Model = typename Converter::OutputModelType;

  const Model& GetOutputModel() const { return GetCvt().GetOutputModel(); }
  Model& GetOutputModel() { return GetCvt().GetOutputModel(); }

  void DeclareAndReportIntSuffix(const SuffixDef<int>& suf,
                                 const std::vector<int>& values) override {
    GetOutputModel().DeclareAndReportIntSuffix(suf, values);
  }


  void DeclareAndReportDblSuffix(const SuffixDef<double>& suf,
    const std::vector<double>& values) override {
    GetOutputModel().DeclareAndReportDblSuffix(suf, values);
  }

  void HandleSolution(int status, fmt::CStringRef msg,
                      const double *x, const double * y, double obj) override {
    GetCvt().HandleSolution(status, msg, x, y, obj);
  }

};

} // namespace mp

#endif // CONVERTER_FLAT_QUERY_H
