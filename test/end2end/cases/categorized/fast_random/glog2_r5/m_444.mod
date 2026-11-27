var x{1..3} integer >= -2 <= 3;

minimize o: x[1] + 6*x[2] + 2*x[3];

s.t. c1: (x[3] == 0) and (x[2] == 1) ==> (x[1] >= 0);
s.t. c2: (if x[1] >= 1 then x[2] else x[2]) < -2;
s.t. c3: max(x[1], x[1]) - min(x[2], x[2]) != 3;
s.t. c4: numberof -1 in (x[1], x[2], x[3]) == 1;
