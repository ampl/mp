var x{1..4} integer >= 0 <= 2;

minimize o: 2*x[1] + 6*x[2] + 5*x[3] + 0*x[4];

s.t. c1: (x[3] != 2) ==> (x[2] != 2);
s.t. c2: (if x[3] > 2 then x[3] else x[1]) != 2;
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4]) < 0;
