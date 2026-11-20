var x{1..3} integer >= 0 <= 2;

maximize o: 3*x[1] + 4*x[2] - 4*x[3];

s.t. c1: exists {i in 1..3} x[i] == 1;
s.t. c2: x[2] == 1;
s.t. c3: x[1] == 2;
