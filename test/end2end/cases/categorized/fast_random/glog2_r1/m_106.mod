var x{1..3} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: (x[2] <= 3) <==> (x[1] + x[1] != 0);
s.t. c2: (if x[2] != 3 then x[3] else x[2]) < 0;
