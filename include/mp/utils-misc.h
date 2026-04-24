#ifndef MP__UTILS_MISC__
#define MP__UTILS_MISC__

namespace mp {

/// Set provided value in destructor
template <class Value>
class RAIIValueSetter {
  Value& where_;
  Value val_;
public:
  /// Construct
  RAIIValueSetter(Value& w, Value v) : where_(w), val_(v) { }
  /// Destruct and set value
  ~RAIIValueSetter() { where_ = std::move(val_); }
};

}

#endif  // MP__UTILS_MISC__
