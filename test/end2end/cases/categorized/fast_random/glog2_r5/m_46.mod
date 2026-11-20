var x{1..5} integer >= -2 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: (x[1] < 1) ==> (x[3] >= 1);
s.t. c2: max(x[3], x[5]) - min(x[5], x[2]) != -2;
s.t. c3: (if x[1] < -2 then x[1] else x[2]) == -1;
s.t. c4: (x[3] > -1) <==> (x[4] + x[4] < 1);
