var x{1..4} integer >= 0 <= 2;

maximize o: -3*x[1] - 3*x[2] + 2*x[3] - 2*x[4];

s.t. c1: ((x[3] <= 0) ==> (x[3] >= 2)) and ((x[3] >= 2) ==> (x[2] >= 0));
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4]) <= 0;
