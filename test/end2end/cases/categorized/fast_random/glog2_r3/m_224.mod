var x{1..5} integer >= -2 <= 2;

minimize o: -2*x[1] + 0*x[2] + 4*x[3] + 6*x[4] + 5*x[5];

s.t. c1: (if x[1] < 2 then x[1] else x[3]) != 0;
