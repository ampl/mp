var x{1..5} integer >= -2 <= 3;

maximize o: 2*x[1] + 3*x[2] - 2*x[3] - 3*x[4] - 5*x[5];

s.t. c1: max(x[4], x[2]) - min(x[3], x[5]) >= -1;
s.t. c2: exists {i in 1..5} x[i] == -2;
s.t. c3: not ((x[2] <= -1) and (x[1] < 2));
