var x{1..6} integer >= 0 <= 1;

maximize o: 2*x[1] - 2*x[2] + 2*x[3] - 4*x[4] - 5*x[5] + 2*x[6];

s.t. c1: exists {i in 1..6} x[i] == 0;
s.t. c2: atleast 3 {i in 1..6} (x[i] >= 1);
