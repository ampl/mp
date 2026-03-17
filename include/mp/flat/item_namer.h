#ifndef ITEM_NAMER_H
#define ITEM_NAMER_H

#include <vector>
#include <string>

namespace mp {

/// Var/con/obj namer for flat models
class ItemNamer {
public:
  /// Default constr
  ItemNamer(std::string nmd) : nm_dfl_(nmd) { }
  /// construct from given names array
  ItemNamer(const std::vector<std::string>& nm, std::string nmd)
      : p_names_given_(&nm), nm_dfl_(nmd) { }
  /// Obtain name[i]
  const char* at(size_t i) {
    if (p_names_given_ && i<p_names_given_->size())
      return (*p_names_given_)[i].c_str();
    if (i>=names_gen_.size())
      names_gen_.resize((size_t)(1.3*i+100));
    if (names_gen_[i].empty())
      names_gen_[i] = nm_dfl_ + std::to_string(i+1) + "_";
    return names_gen_[i].c_str();
  }
private:
  const std::vector<std::string>* p_names_given_ {};
  std::vector<std::string> names_gen_;
  std::string nm_dfl_;
};

}  // namespace mp

#endif // ITEM_NAMER_H
