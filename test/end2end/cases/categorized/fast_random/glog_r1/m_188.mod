var x{1..3} integer >= -2 <= 0;

minimize o: -3*x[1] + 0*x[2] + 3*x[3];

s.t. c1: (if x[1] < -1 then x[1] else x[3]) > -1;
