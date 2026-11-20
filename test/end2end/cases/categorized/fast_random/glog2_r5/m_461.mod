var x{1..4} integer >= -2 <= 2;

maximize o: -5*x[1] + 3*x[2] - 3*x[3] - 5*x[4];

s.t. c1: not ((x[3] == -3) and (x[3] >= 2));
s.t. c2: (x[2] == 3) and (x[3] >= 3) ==> (x[1] >= 3);
s.t. c3: (x[4] != 2) or (x[4] > 2) or (x[3] != 3);
