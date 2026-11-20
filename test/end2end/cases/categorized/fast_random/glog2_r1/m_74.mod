var x{1..5} integer >= -2 <= 2;

minimize o: -x[1] + 2*x[2] - 3*x[3] - 3*x[4] + 6*x[5];

s.t. c1: (x[4] < -2) or (x[3] <= -2) or (x[3] < -1);
s.t. c2: (x[5] > -1) ==> (x[2] >= 1) else (x[2] != -2);
s.t. c3: (if x[4] < -2 then x[1] else x[2]) > 2;
