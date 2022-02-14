#ifndef CONVERTERBASE_H
#define CONVERTERBASE_H

#include "mp/env.h"

namespace mp {

/// BasicConverter<> is an interface for a class that takes
/// a ProblemBuilder and translates it into the solver API.
/// The ProblemBuilder can be used by NLHandler for NL input.
template <class ProblemBuilder>
class BasicConverter : public EnvKeeper {
public:
  /// Constructor
  BasicConverter(Env& e) : EnvKeeper(e) { }
  /// Destructor
  virtual ~BasicConverter() = default;

  /// The ProblemBuilder a.k.a. Model type
  using ModelType = ProblemBuilder;

  /// Access the model instance, const
  virtual const ModelType& GetModel() const = 0;
  /// Access the model instance
  virtual ModelType& GetModel() = 0;

  /// Init solver options
  virtual void InitOptions() { }

  /// Convert the model into the solver API
  virtual void ConvertModel() = 0;
};

} // namespace mp

#endif // CONVERTERBASE_H
