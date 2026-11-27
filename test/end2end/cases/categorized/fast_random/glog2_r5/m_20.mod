var x{1..3} integer >= -2 <= 0;

minimize o: -2*x[1] - 2*x[2] - 5*x[3];

s.t. c1: (if x[2] <= -1 then x[3] else x[1]) >= 0;
