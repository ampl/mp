var x{1..5} integer >= 0 <= 4;

maximize o: -2*x[1] + 4*x[2] - 4*x[3] + 0*x[4] - 3*x[5];

s.t. c1: max(x[1], x[4]) - min(x[5], x[3]) >= 0;
s.t. c2: (if x[2] == 4 then x[5] else x[4]) < -1;
