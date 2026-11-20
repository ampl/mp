var x{1..3} integer >= -2 <= 1;

minimize o: 6*x[1] - 5*x[2] + 2*x[3];

s.t. c1: (x[2] == -1) <==> (x[3] + x[3] != -2);
s.t. c2: (if x[3] <= -2 then x[3] else x[3]) != 0;
