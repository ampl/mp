var x{1..3} integer >= 0 <= 4;

maximize o: 5*x[1] + 2*x[2] - x[3];

s.t. c1: max(x[2], x[3]) - min(x[1], x[3]) < 2;
s.t. c2: exists {i in 1..3} x[i] == 4;
s.t. c3: x[2] == 4;
