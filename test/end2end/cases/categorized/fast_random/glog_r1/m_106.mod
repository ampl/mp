var x{1..5} integer >= 0 <= 2;

maximize o: -5*x[1] - 2*x[2] - 2*x[3] - 3*x[4] + 4*x[5];

s.t. c1: atleast 1 {i in 1..5} (x[i] >= 2);
s.t. c2: exists {i in 1..5} x[i] == 2;
