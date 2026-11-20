var x{1..4} integer >= 0 <= 3;

maximize o: 2*x[1] + 3*x[2] + 3*x[3] - 4*x[4];

s.t. c1: ((x[4] > 2) ==> (x[1] != 1)) and ((x[1] != 1) ==> (x[2] >= 2));
s.t. c2: not ((x[2] >= 3) and (x[2] > 3));
s.t. c3: max(x[2], x[2]) - min(x[1], x[2]) < 3;
