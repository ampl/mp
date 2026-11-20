var x{1..3} integer >= 0 <= 3;

minimize o: -3000*x[1] + 17*x[2] + 17*x[3];

s.t. c1: (if x[2] >= 0 then x[2] else x[2]) < 2;
s.t. c2: (x[3] >= 1) <==> (x[1] + x[1] != 1);
