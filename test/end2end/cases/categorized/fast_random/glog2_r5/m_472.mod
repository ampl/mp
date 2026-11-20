var x{1..3} integer >= -2 <= 1;

minimize o: -2*x[1] + 3*x[2] - 4*x[3];

s.t. c1: abs(x[2] - x[3]) != 0;
s.t. c2: (if x[1] == 0 then x[1] else x[1]) > 0;
s.t. c3: (x[2] >= -2) and (x[3] == 1) ==> (x[2] >= -2);
s.t. c4: (x[3] <= -1) or (x[3] >= 1) or (x[3] >= 1);
