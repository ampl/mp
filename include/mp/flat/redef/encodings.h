#ifndef ENCODINGS_H
#define ENCODINGS_H

#include <vector>
#include <memory>

namespace mp {

class ZZI_Encoding;

/// Obtain an object of ZII_Encoding
using P_ZZI_Encoding = std::shared_ptr<ZZI_Encoding>;
P_ZZI_Encoding MakeZZIEncoding();

/// Obtain column \a k of the encoding matrix C^r
/// prepended by 0 and postpended by line v1-1
/// if v1-1 == 2^r
std::vector<double>
GetExtendedColumn(ZZI_Encoding& ,
                  int r, int k, int v0, int v1);

}  // namespace mp

#endif // ENCODINGS_H
