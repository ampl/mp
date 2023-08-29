#ifndef MP_UTILS_MATH_H
#define MP_UTILS_MATH_H

#include <cmath>

namespace mp {

/// https://stackoverflow.com/questions/13094224/a-c-routine-to-round-a-float-to-n-significant-digits
template <class F>
F round_to_digits(F value, int digits) {
    if (value == 0.0) // otherwise it will return 'nan'
        return 0.0;   // due to the log10() of zero

    F factor = std::pow(10.0,
                        digits
                        - std::ceil(std::log10(std::fabs(value))));
    return std::round(value * factor) / factor;
}

}  // namespace mp

#endif // MP_UTILS_MATH_H
