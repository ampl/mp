var x{1..4} integer >= 0 <= 2;

maximize o: -5*x[1] + 4*x[2] + 2*x[3] + x[4];

s.t. c1: (if x[2] != 1 then x[2] else x[1]) == 1;
