var x{1..3} integer >= 0 <= 3;

maximize o: 5*x[1] + 4*x[2] - 4*x[3];

s.t. c1: abs(x[1] - x[1]) != 0;
s.t. c2: (if x[2] <= 0 then x[2] else x[1]) == 3;
