var x{1..4} integer >= 0 <= 3;

minimize o: -2*x[1] - 5*x[2] - 2*x[3] - x[4];

s.t. c1: not (x[2] < 3);
s.t. c2: (if x[2] <= 3 then x[2] else x[1]) < 2;
