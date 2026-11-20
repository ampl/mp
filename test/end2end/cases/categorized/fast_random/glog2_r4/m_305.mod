var x{1..3} integer >= 0 <= 2;

minimize o: 0*x[1] + 5*x[2] - 3*x[3];

s.t. c1: (x[1] == 0) ==> (x[2] >= 1);
s.t. c2: (x[3] != 1) and (x[3] == 0) ==> (x[2] >= 1);
s.t. c3: max(x[1], x[1]) - min(x[3], x[2]) == 2;
s.t. c4: abs(x[1] - x[1]) == 1;
s.t. c5: numberof 0 in (x[1], x[2], x[3]) >= 3;
