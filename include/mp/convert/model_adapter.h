#ifndef MODEL_ADAPTER_H
#define MODEL_ADAPTER_H

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


  typename Model::IntSuffixHandler
  AddIntSuffix(fmt::StringRef name, int kind, int=0) {
    return Model::AddIntSuffix(name, kind);
  }

};

#endif // MODEL_ADAPTER_H
