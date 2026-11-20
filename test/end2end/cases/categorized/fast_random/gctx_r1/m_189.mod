var x{1..5} integer >= 0 <= 2;

maximize o: 4*x[1] - x[2] - 5*x[3] + 2*x[4] + 3*x[5];

s.t. c1: (if abs(x[5] - x[3]) < 0 then x[4] else x[3]) <= 0;
s.t. c2: not ((abs(x[5] - x[3]) < 0) and (x[1] != 2));
