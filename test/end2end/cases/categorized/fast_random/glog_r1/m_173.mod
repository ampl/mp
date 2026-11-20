var x{1..6} integer >= -2 <= 1;

minimize o: 3*x[1] + 0*x[2] - 2*x[3] + 6*x[4] - 5*x[5] + 2*x[6];

s.t. c1: (if x[1] <= -1 then x[5] else x[5]) < 0;
s.t. c2: (x[5] == -2) <==> (x[1] + x[4] == 1);
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 2;
