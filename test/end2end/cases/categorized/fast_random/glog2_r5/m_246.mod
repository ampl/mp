var x{1..4} integer >= -2 <= 0;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: ((x[1] > -2) ==> (x[4] < 0)) and ((x[4] < 0) ==> (x[4] >= -2));
