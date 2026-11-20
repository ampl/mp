var x{1..6} integer >= -2 <= 1;

minimize o: -5*x[1] - x[2] + 4*x[3] + 3*x[4] + 4*x[5] + 4*x[6];

s.t. c1: (x[5] <= -1) or (x[5] != 0) or (x[1] <= -2);
s.t. c2: (x[4] <= 1) ==> (x[3] != -1) else (x[3] == 1);
s.t. c3: (x[6] <= -2) or (x[2] < 0) or (x[3] < 1);
s.t. c4: abs(x[2] - x[2]) <= 0;
s.t. c5: numberof -1 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 1;
s.t. c6: abs(x[2] - x[6]) < 3;
