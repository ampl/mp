#ifndef MP2MIP_H
#define MP2MIP_H

#include "mp/convert/converter_flat.h"

namespace mp {

/// MPToMIPConverter: one of the converters requiring a "minimal" output interface
template <class Impl, class Backend,
          class Model = BasicModel<std::allocator<char> > >
class MPToMIPConverter
    : public BasicMPFlatConverter<Impl, Backend, Model>
{

};

} // namespace mp

#endif // MP2MIP_H
