var x{1..3} integer >= 0 <= 2;

maximize o: 5*x[1] - 4*x[2] + 5*x[3];

s.t. c1: (if x[1] != 1 then x[1] else x[1]) > 1;
