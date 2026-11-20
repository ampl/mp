var x{1..3} integer >= 0 <= 4;

maximize o: 3*x[1] - x[2] + 0*x[3];

s.t. c1: (x[1] > 2) or (x[1] < 4) or (x[3] >= 0);
s.t. c2: not ((x[3] != 1) and (x[1] <= 3));
s.t. c3: max(x[1], x[2]) - min(x[1], x[3]) < 4;
s.t. c4: (x[3] >= 3) ==> (x[3] >= 3);
s.t. c5: max(x[2], x[3]) - min(x[3], x[1]) != 0;
s.t. c6: (x[2] < 2) and (x[3] > 2) ==> (x[3] >= 2);
