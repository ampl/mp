var x{1..5} integer >= 0 <= 4;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: (if x[1] == 4 then x[4] else x[3]) != 2;
s.t. c2: abs(x[3] - x[1]) < 1;
