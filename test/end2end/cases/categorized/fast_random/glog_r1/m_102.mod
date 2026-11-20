var x{1..3} integer >= -2 <= 2;

minimize o: -4*x[1] - 4*x[2] - 4*x[3];

s.t. c1: numberof 0 in (x[1], x[2], x[3]) != 2;
s.t. c2: (x[2] != 2) and (x[1] >= -2) ==> (x[3] >= 2);
s.t. c3: (if x[1] >= 0 then x[2] else x[1]) < -1;
s.t. c4: (x[2] <= 1) <==> (x[1] + x[1] >= 0);
