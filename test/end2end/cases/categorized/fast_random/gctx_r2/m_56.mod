var x{1..4} integer >= -1 <= 1;

minimize o: 2*x[1] - 3*x[2] + 4*x[3] + 6*x[4];

s.t. c1: (min(x[2], x[4]) <= -1) <==> (x[4] >= 1);
s.t. c2: not (min(x[4], x[4]) < 1);
s.t. c3: (x[4] >= 1) ==> ((x[4] >= 1) or (min(x[2], x[4]) <= -1));
s.t. c4: not (x[4] >= 1);
