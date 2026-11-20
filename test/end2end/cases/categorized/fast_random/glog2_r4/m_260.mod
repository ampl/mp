var x{1..3} integer >= -2 <= 1;

minimize o: -5*x[1] + 3*x[2] + 2*x[3];

s.t. c1: (if x[2] <= -2 then x[2] else x[2]) >= -1;
