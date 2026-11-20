var x{1..3} integer >= -2 <= 2;

minimize o: -5*x[1] - 2*x[2] - 4*x[3];

s.t. c1: (if x[3] != 1 then x[2] else x[1]) == -1;
