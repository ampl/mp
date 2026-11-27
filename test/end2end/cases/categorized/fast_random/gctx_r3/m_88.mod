var x{1..4} integer >= 0 <= 2;

maximize o: 6*x[1] - 5*x[2] - x[3] + 4*x[4];

s.t. c1: (x[3] < 0) ==> (max(x[1], x[1]) <= 1);
s.t. c2: (x[3] < 0) or (x[3] < 0);
s.t. c3: (max(x[1], x[1]) <= 1) ==> ((max(x[1], x[1]) <= 1) or (x[3] < 0));
