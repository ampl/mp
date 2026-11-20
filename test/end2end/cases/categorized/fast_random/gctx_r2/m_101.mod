var x{1..5} integer >= 0 <= 3;

minimize o: 5*x[1] + 0*x[2] - x[3] + 0*x[4] + 6*x[5];

s.t. c1: (if x[4] < 1 then x[1] else 0) + (if x[3] == 0 then x[5] else 0) < 1;
s.t. c2: (if x[3] == 0 then x[5] else 0) + (if x[3] == 0 then x[1] else 0) >= 3;
