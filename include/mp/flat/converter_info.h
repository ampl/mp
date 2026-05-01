#ifndef CONVERTER_INFO_H
#define CONVERTER_INFO_H

namespace mp {

/// Information from a model converter
class ConverterInfo {
public:
  /// Destruct
  virtual ~ConverterInfo() { }

  /// Chosen value of the maximal refcount
  /// for an algebraic subexpression
  virtual int RefCountMaxAlgebraic() const = 0;
};

}  // namespace mp

#endif // CONVERTER_INFO_H
