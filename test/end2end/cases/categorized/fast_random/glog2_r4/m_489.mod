var x{1..4} integer >= 0 <= 2;

minimize o: -5*x[1] + 3*x[2] + 4*x[3] + x[4];

s.t. c1: (if x[1] == 0 then x[2] else x[3]) == 0;
s.t. c2: (if x[1] == 0 then x[2] else x[3]) == 0;
