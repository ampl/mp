var x{1..5} integer >= 0 <= 4;

minimize o: -5*x[1] - 2*x[2] + 2*x[3] + 0*x[4] - 2*x[5];

s.t. c1: (x[2] < 2) or (x[1] > 2) or (x[5] <= 1);
s.t. c2: ((if x[4] <= 4 then 1 else 0) + (if x[4] >= 4 then 1 else 0)) >= 1;
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) >= 3;
s.t. c4: x[2] == 0;
s.t. c5: x[3] == 3;
