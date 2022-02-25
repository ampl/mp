#ifndef REDEF_BASE_H
#define REDEF_BASE_H

#include "mp/format.h"
#include "mp/presolve.h"

namespace mp {

/// A base class for specific Item Converters, such as
/// individual constraint converters.
/// @param ModelConverter: should have methods
/// GetPresolver() and GetConstraintKeeper(Constraint*).
template <class ModelConverter>
class BasicItemConverter {
public:
  /// Constructor
  BasicItemConverter(ModelConverter& mc) : mdl_cvt_(mc) { }

  /// Access ModelConverter
  ModelConverter& GetMC() { return mdl_cvt_; }
  /// Access Presolver
  pre::Presolver& GetPre() { return GetMC().GetPresolver(); }
  /// Access specific item's ValueNode
  template <class Item>
  pre::ValueNode& GetVN(Item* pc)
  { return GetMC().GetConstraintKeeper(pc).GetValueNode(); }

private:
  ModelConverter& mdl_cvt_;
};


/// A functional constraint converter
/// @param Impl: the final converter
/// @param ModelConverter
template <class Impl, class ModelConverter>
class BasicFuncConstrCvt :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  BasicFuncConstrCvt(ModelConverter& mc) : Base(mc) { }

  /// Generic Convert(), distinguishes context & result variable
  ///
  /// Responsible for adding presolve bridges, if any
  /// @param item: the item to be converted
  /// @param i: item index, used to create a presolve bridge
  ///
  /// The Impl can reimplement
  template <class ItemType>
  void Convert(const ItemType& item, int i) {
    auto ctx = item.GetContext();
    auto rv = item.GetResultVar();
    auto bnd00 = item.GetAprioriBounds();   // TODO use current bounds?
    if ( ctx.HasNegative() &&
         GetMC().lb(rv) < bnd00.second ) {  // Need the negative direction
      MPD( ConvertCtxNeg(item, i) );
    }
    if ( ctx.HasPositive() &&
         GetMC().ub(rv) > bnd00.first ) {   // Need the positive direction
      MPD( ConvertCtxPos(item, i) );
    }
    // TODO infeasible cases?
  }

  /// Convert in negative context
  template <class ItemType>
  void ConvertCtxNeg(const ItemType& item, int ) {
    MP_RAISE( fmt::format(
                "Conversion of '{}' in negative context "
                "not implemented", item.GetName() ) );
  }
  /// Convert in positive context
  template <class ItemType>
  void ConvertCtxPos(const ItemType& item, int ) {
    MP_RAISE( fmt::format(
                "Conversion of '{}' in positive context "
                "not implemented", item.GetName() ) );
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// To be used by descendants of BasiccOnverter
#define GET_CONSTRAINT_VALUE_NODE(con_type) \
  this->GetMC().GetValueNode((con_type*)nullptr)

/// In the ModelConverter: to use a specific item_cvt_type<>
/// Assumes Impl is the final ModelConverter type
#define INSTALL_ITEM_CONVERTER(item_cvt_type) \
  item_cvt_type<Impl> item_cvt__ ## item_cvt_type ## _ \
   { *static_cast<Impl*>(this) }; \
  void Convert(const typename \
      item_cvt_type<Impl>::ItemType& con, int i) { \
    item_cvt__ ## item_cvt_type ## _ . Convert(con, i); \
  }


} // namespace mp

#endif // REDEF_BASE_H
