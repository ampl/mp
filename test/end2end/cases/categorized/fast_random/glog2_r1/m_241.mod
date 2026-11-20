var x{1..5} integer >= 0 <= 4;

minimize o: -3*x[1] + 2*x[2] + 5*x[3] + 5*x[4] + 5*x[5];

s.t. c1: atleast 3 {i in 1..5} (x[i] >= 5);
s.t. c2: exists {i in 1..5} x[i] == 5;
