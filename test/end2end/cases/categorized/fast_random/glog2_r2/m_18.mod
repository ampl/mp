var x{1..7} integer >= 0 <= 1;

maximize o: x[1] - 4*x[2] - 2*x[3] + 6*x[4] - 2*x[5] - 3*x[6] - 5*x[7];

s.t. c1: (if x[1] == 0 then x[2] else x[1]) > 0;
