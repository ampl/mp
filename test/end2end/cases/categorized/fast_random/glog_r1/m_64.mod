var x{1..4} integer >= 0 <= 3;

maximize o: x[1] + x[2] + 4*x[3] + x[4];

s.t. c1: (x[4] > 3) and (x[2] == 3) ==> (x[4] >= 3);
s.t. c2: (x[1] == 1) and (x[1] > 2) ==> (x[3] >= 1);
