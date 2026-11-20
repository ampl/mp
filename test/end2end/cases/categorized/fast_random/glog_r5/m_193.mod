var x{1..3} integer >= 0 <= 4;

maximize o: 2*x[1] + x[2] - 3*x[3];

s.t. c1: (if x[1] >= 3 then x[1] else x[1]) == 1;
