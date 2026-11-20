var x{1..4} integer >= 0 <= 2;

minimize o: -x[1] + 2*x[2] + 0*x[3] + 6*x[4];

s.t. c1: (x[4] > 2) ==> ((x[4] > 2) or (x[4] > 2));
s.t. c2: (x[4] + x[1] <= 2) and ((x[4] > 2) or not (x[4] + x[1] <= 2));
