var x{1..4} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (if x[2] > 2 then x[3] else x[4]) < 4;
