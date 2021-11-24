/*
 A mathematical optimization solver.

 Copyright (C) 2021 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */

#include "mp/flat/std_constr.h"

mp::PLPoints::PLPoints(const PLSlopes &pls) {
  constexpr auto eps = 1.0;           // works?
  const auto& bp = pls.GetBP();
  const auto& sl = pls.GetSlopes();
  const auto nsl = sl.size();
  const auto X0 = pls.GetX0(), Y0 = pls.GetY0();
  x_.resize(nsl+1);
  y_.resize(nsl+1);
  /// Copy and add dummy breakpoints on both ends
  std::copy(bp.begin(), bp.end(), x_.begin()+1);
  x_[0] = x_[1] - eps;
  x_[nsl] = x_[nsl-1] + eps;
  y_[0] = 0.0;                        // initialize leftmost point
  /// Lift the line by this
  double deltaH {};
  if (x_[0] > X0)                     // if left x > reference point
    deltaH = sl[0] * (x_[0]-X0) + Y0;
  for (size_t i = 0; i < nsl; ++i) {
    assert( x_[i+1] > x_[i] );
    y_[i+1] = y_[i] + sl[i] * (x_[i+1]-x_[i]);
    if (x_[i]<=X0 && (x_[i+1]>=X0 || i==nsl-1))
      deltaH = Y0 - (y_[i] + sl[i]*(-x_[i]-X0));
  }
  for (size_t i = 0; i <= nsl; ++i) {
    y_[i] += deltaH;
  }
}
