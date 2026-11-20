var x{1..4} integer >= 0 <= 4;

maximize o: -2*x[1] - 4*x[2] + 2*x[3] - 3*x[4];

s.t. c1: (x[1] < 3) ==> (x[2] > 2);
s.t. c2: x[3] == 4;
