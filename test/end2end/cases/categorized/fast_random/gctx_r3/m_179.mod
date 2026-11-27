var x{1..3} integer >= 0 <= 3;

minimize o: 5*x[1] - 3*x[2] + 2*x[3];

s.t. c1: (x[3] >= 3) ==> (x[3] >= 3);
s.t. c2: not (x[3] >= 3);
s.t. c3: (if x[3] >= 3 then x[2] else x[3]) <= 0;
s.t. c4: (if x[1] > 2 then x[1] else x[2]) == 2;
