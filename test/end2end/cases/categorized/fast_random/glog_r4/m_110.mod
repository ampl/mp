var x{1..3} integer >= -2 <= 2;

maximize o: -x[1] - 3*x[2] - x[3];

s.t. c1: (if x[2] <= 1 then x[3] else x[1]) <= 0;
