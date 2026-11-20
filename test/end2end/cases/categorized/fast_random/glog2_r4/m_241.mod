var x{1..4} integer >= 0 <= 3;

maximize o: 4*x[1] + 0*x[2] - 4*x[3] - 4*x[4];

s.t. c1: (if x[4] > 3 then x[3] else x[2]) != -1;
s.t. c2: abs(x[2] - x[1]) >= 3;
s.t. c3: abs(x[2] - x[2]) > 2;
