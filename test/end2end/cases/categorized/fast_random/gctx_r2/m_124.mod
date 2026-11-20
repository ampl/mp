var x{1..5} integer >= 0 <= 2;

maximize o: -5*x[1] - 2*x[2] + 0*x[3] - 3*x[4] - x[5];

s.t. c1: not (abs(x[3] - x[4]) < 0);
s.t. c2: (min(x[5], x[3]) > 1) ==> (abs(x[3] - x[4]) < 0);
s.t. c3: not ((x[4] == x[3]) and (min(x[5], x[3]) > 1));
s.t. c4: not (x[4] == x[3]);
