var x{1..3} integer >= -2 <= 2;

maximize o: -4*x[1] + 4*x[2] - 3*x[3];

s.t. c1: (x[2] > -1) ==> (x[2] < -2) else (x[2] >= -2);
s.t. c2: numberof 1 in (x[1], x[2], x[3]) <= 1;
