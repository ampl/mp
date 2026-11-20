var x{1..3} integer >= -2 <= 0;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: ((x[1] > 0) ==> (x[1] < -1)) and ((x[1] < -1) ==> (x[2] >= 0));
s.t. c2: abs(x[3] - x[3]) == 0;
s.t. c3: (if x[2] >= -1 then x[3] else x[3]) != -1;
s.t. c4: max(x[2], x[2]) - min(x[3], x[2]) == 0;
