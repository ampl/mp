var x{1..3} integer >= -3 <= 1;

maximize o: 2*x[1] + 0*x[2] + 4*x[3];

s.t. c1: (x[3] >= 0) ==> (x[2] > -1) else (x[3] != -3);
s.t. c2: (x[1] <= 1) or (x[2] > -2) or (x[3] == -1);
s.t. c3: (if x[2] >= -3 then x[1] else x[1]) != -3;
