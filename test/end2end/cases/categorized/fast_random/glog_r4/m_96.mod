var x{1..3} integer >= 0 <= 3;

minimize o: 2*x[1] + x[2] - 2*x[3];

s.t. c1: (if x[2] <= 2 then x[2] else x[3]) == 0;
