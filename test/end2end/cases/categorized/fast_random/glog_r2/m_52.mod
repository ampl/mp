var x{1..6} integer >= -2 <= 0;

minimize o: 5*x[1] + 6*x[2] - 4*x[3] - x[4] - 4*x[5] + 0*x[6];

s.t. c1: (x[1] >= 0) or (x[4] >= 0) or (x[1] != -2);
s.t. c2: (if x[3] >= 0 then x[3] else x[1]) == -2;
s.t. c3: (if x[5] <= -1 then x[4] else x[4]) < -1;
s.t. c4: (x[5] <= -1) and (x[2] == -1) ==> (x[6] >= -1);
