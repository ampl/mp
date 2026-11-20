var x{1..3} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: numberof 2 in (x[1], x[2], x[3]) > 3;
s.t. c2: (x[3] >= 3) <==> (x[2] + x[1] == 0);
s.t. c3: (x[1] != 3) and (x[3] == 0) ==> (x[1] >= 3);
s.t. c4: abs(x[1] - x[2]) != 2;
