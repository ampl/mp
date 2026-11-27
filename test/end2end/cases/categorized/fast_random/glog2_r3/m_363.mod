var x{1..6} integer >= 0 <= 4;

minimize o: -x[1] + 6*x[2] + 6*x[3] + 4*x[4] + x[5] - 2*x[6];

s.t. c1: (x[2] >= 4) and (x[1] > 1) ==> (x[2] >= 4);
s.t. c2: max(x[5], x[5]) - min(x[6], x[5]) > 1;
s.t. c3: (x[1] < 2) and (x[5] <= 3) ==> (x[4] >= 2);
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 0;
s.t. c5: x[6] == 1;
