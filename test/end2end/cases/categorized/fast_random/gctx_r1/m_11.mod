var x{1..5} integer >= 0 <= 2;

maximize o: -x[1] + 2*x[2] + 2*x[3] + 0*x[4] + 2*x[5];

s.t. c1: (min(x[5], x[3]) == 2) or (x[3] == 2);
s.t. c2: not (min(x[4], x[5]) > 0);
