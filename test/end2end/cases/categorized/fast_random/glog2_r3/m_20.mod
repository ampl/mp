var x{1..4} integer >= 0 <= 3;

maximize o: 3*x[1] + 0*x[2] + 0*x[3] - x[4];

s.t. c1: (x[2] != 2) or (x[2] == 1) or (x[1] > 2);
s.t. c2: ((x[3] > 0) ==> (x[4] == 3)) and ((x[4] == 3) ==> (x[3] >= 0));
s.t. c3: not ((x[2] <= 3) and (x[2] >= 0));
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4]) >= 2;
