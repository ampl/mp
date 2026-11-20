var x{1..3} integer >= -2 <= 0;

maximize o: -4*x[1] + 3*x[2] + 2*x[3];

s.t. c1: alldiff(x[1], x[3]);
s.t. c2: max(x[3], x[1]) - min(x[3], x[2]) >= -2;
s.t. c3: exists {i in 1..3} x[i] == -2;
s.t. c4: alldiff(x[1], x[3]);
s.t. c5: exists {i in 1..3} x[i] == -2;
