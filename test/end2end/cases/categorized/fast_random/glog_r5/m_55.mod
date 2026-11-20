var x{1..4} integer >= 0 <= 3;

maximize o: -4*x[1] + 5*x[2] - 3*x[3] + x[4];

s.t. c1: ((x[3] != 0) ==> (x[1] < 0)) and ((x[1] < 0) ==> (x[1] >= 0));
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) >= 0;
s.t. c3: not ((x[3] != 0) and (x[1] != 1));
s.t. c4: (x[2] != 3) ==> (x[1] <= 1);
