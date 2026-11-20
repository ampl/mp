var x{1..5} integer >= 0 <= 3;

minimize o: 6*x[1] - 4*x[2] - 5*x[3] + 2*x[4] - 4*x[5];

s.t. c1: (if x[4] >= -1 then x[1] else x[5]) <= 4;
s.t. c2: abs(x[3] - x[2]) == 3;
s.t. c3: (x[1] >= -1) and (x[1] == 4) ==> (x[5] >= -1);
