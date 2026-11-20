var x{1..3} integer >= 0 <= 2;

maximize o: -4*x[1] + 5*x[2] - 3*x[3];

s.t. c1: atleast 1 {i in 1..3} (x[i] >= 1);
s.t. c2: not ((x[3] != 2) and (x[2] > 1));
