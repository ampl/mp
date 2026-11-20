var x{1..3} integer >= -2 <= 1;

minimize o: 0*x[1] + 4*x[2] - 4*x[3];

s.t. c1: (x[1] < 0) or (x[2] < 0) or (x[2] <= 1);
s.t. c2: (x[1] >= 0) and (x[3] < -2) ==> (x[2] >= 0);
s.t. c3: max(x[2], x[3]) - min(x[1], x[1]) > 0;
s.t. c4: (if x[2] != 1 then x[1] else x[1]) >= 0;
