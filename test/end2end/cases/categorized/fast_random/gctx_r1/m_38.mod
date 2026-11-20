var x{1..5} integer >= -1 <= 2;

maximize o: 4*x[1] - 2*x[2] + 2*x[3] + 0*x[4] + 6*x[5];

s.t. c1: not (min(x[2], x[5]) <= -1);
s.t. c2: (if min(x[2], x[5]) <= -1 then x[2] else x[4]) > 2;
