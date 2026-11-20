var x{1..6} integer >= 0 <= 3;

maximize o: 2*x[1] + 4*x[2] + 4*x[3] - 4*x[4] + 2*x[5] + x[6];

s.t. c1: max(x[2], x[3]) - min(x[6], x[1]) >= 1;
s.t. c2: ((x[3] != 3) ==> (x[6] >= 0)) and ((x[6] >= 0) ==> (x[2] >= 3));
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 2;
