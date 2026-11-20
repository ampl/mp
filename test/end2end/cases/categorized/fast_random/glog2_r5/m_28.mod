var x{1..7} integer >= 0 <= 1;

minimize o: -5*x[1] - x[2] + 6*x[3] - 2*x[4] - 4*x[5] + x[6] + 2*x[7];

s.t. c1: (x[5] <= 0) ==> (x[5] <= 1);
s.t. c2: (if x[3] == 1 then x[5] else x[1]) == 0;
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6], x[7]) >= 3;
