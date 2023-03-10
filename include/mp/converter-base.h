#ifndef CONVERTERBASE_H
#define CONVERTERBASE_H

#include "mp/env.h"
#include "mp/ampls-ccallbacks.h"

namespace mp {

/// BasicConverter<> is an interface for a class that takes
/// a Model and translates it into the solver API.
/// Model can be a ProblemBuilder used by NLHandler for NL input,
/// as done by ModelManagerWithProblemBuilder.
template <class Model>
class BasicConverter : public EnvKeeper {
public:
  /// Constructor
  BasicConverter(Env& e) : EnvKeeper(e) { }
  /// Destructor
  virtual ~BasicConverter() = default;

  /// The ProblemBuilder a.k.a. Model type
  using ModelType = Model;

  /// Access the model instance, const
  virtual const ModelType& GetModel() const = 0;
  /// Access the model instance
  virtual ModelType& GetModel() = 0;

  /// Init solver options
  virtual void InitOptions() { }

  /// Convert the model into the solver API
  virtual void ConvertModel() = 0;

  /// Fill model traits
  virtual void FillModelTraits(AMPLS_ModelTraits& ) = 0;
};

} // namespace mp

#endif // CONVERTERBASE_H
