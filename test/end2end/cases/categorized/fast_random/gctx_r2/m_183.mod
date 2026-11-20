var x{1..4} integer >= -1 <= 2;

minimize o: -4*x[1] - 5*x[2] + x[3] + 2*x[4];

s.t. c1: (x[1] >= 0) ==> (x[1] >= 0);
s.t. c2: (if x[1] >= 0 then x[1] else x[2]) > 2;
s.t. c3: (x[1] == -1) or (x[1] == -1);
s.t. c4: not (x[1] == -1);
