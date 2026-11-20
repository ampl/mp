var x{1..5} integer >= 0 <= 2;

minimize o: -3*x[1] - x[2] + x[3] + 0*x[4] + 0*x[5];

s.t. c1: ((if x[3] != 2 then 1 else 0) + (if x[5] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[3] == 0) or (x[3] >= 1) or (x[5] != 1);
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) < 0;
