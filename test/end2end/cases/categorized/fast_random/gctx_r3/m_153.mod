var x{1..4} integer >= 0 <= 2;

maximize o: 5*x[1] - 2*x[2] + 4*x[3] + 4*x[4];

s.t. c1: not (min(x[4], x[4]) == 0);
s.t. c2: (x[2] < x[3]) <==> (x[2] < x[3]);
s.t. c3: not (abs(x[1] - x[3]) != 1);
s.t. c4: not (x[2] < x[3]);
