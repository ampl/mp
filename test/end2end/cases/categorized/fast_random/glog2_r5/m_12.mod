var x{1..5} integer >= 0 <= 4;

maximize o: x[1] - 4*x[2] + 5*x[3] - x[4] + 6*x[5];

s.t. c1: atleast 3 {i in 1..5} (x[i] >= 1);
s.t. c2: atmost 2 {i in 1..5} (x[i] >= 1);
s.t. c3: alldiff(x[1], x[2], x[3], x[4]);
