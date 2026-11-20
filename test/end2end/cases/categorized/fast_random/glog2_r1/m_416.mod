var x{1..6} integer >= 0 <= 4;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: (x[3] >= 4) and (x[4] < 0) ==> (x[3] >= 4);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 1;
s.t. c3: max(x[6], x[1]) - min(x[5], x[5]) < 0;
