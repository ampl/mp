var x{1..3} integer >= 0 <= 2;

minimize o: 6*x[1] + x[2] + 6*x[3];

s.t. c1: (if x[3] == 1 then x[2] else x[3]) < 2;
