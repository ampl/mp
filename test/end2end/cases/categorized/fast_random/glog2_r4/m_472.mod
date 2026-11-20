var x{1..6} integer >= -2 <= 0;

minimize o: -4*x[1] + x[2] - 4*x[3] + 5*x[4] + 6*x[5] + 2*x[6];

s.t. c1: (x[5] >= 0) ==> (x[2] > -2);
s.t. c2: ((if x[1] <= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 0;
s.t. c4: (if x[6] < -1 then x[3] else x[6]) == 0;
