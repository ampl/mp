var x{1..5} integer >= 0 <= 3;

minimize o: x[1] + 0*x[2] - 3*x[3] + 5*x[4] + 2*x[5];

s.t. c1: (x[1] > -1) and (x[4] <= 4) ==> (x[3] >= -1);
s.t. c2: abs(x[4] - x[4]) >= 0;
s.t. c3: max(x[1], x[4]) - min(x[2], x[2]) < 4;
s.t. c4: (if x[1] < 0 then x[4] else x[4]) != 0;
