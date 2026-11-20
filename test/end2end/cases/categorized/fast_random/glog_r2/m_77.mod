var x{1..5} integer >= 0 <= 4;

maximize o: 4*x[1] + 0*x[2] + 4*x[3] + 4*x[4] - x[5];

s.t. c1: max(x[2], x[2]) - min(x[5], x[1]) > 2;
s.t. c2: exists {i in 1..5} x[i] == 4;
s.t. c3: ((x[5] == 4) ==> (x[1] != 1)) and ((x[1] != 1) ==> (x[1] >= 4));
