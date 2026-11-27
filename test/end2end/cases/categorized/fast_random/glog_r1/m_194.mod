var x{1..3} integer >= -2 <= 1;

minimize o: 0*x[1] - 5*x[2] + 2*x[3];

s.t. c1: (if x[1] < -2 then x[2] else x[3]) > -1;
s.t. c2: (x[2] == -2) <==> (x[2] + x[1] >= -1);
