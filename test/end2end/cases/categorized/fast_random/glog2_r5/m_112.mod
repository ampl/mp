var x{1..6} integer >= -2 <= 2;

minimize o: 2*x[1] - 4*x[2] + 5*x[3] + 5*x[4] + 2*x[5] + 5*x[6];

s.t. c1: ((if x[4] >= -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c2: not ((x[5] > 0) and (x[5] <= 2));
s.t. c3: abs(x[5] - x[5]) == 0;
s.t. c4: (x[2] != 1) ==> (x[2] < -1);
