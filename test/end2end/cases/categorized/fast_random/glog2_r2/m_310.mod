var x{1..5} integer >= -2 <= 2;

maximize o: 5*x[1] - 5*x[2] - x[3] - 4*x[4] - 5*x[5];

s.t. c1: (if x[1] == -2 then x[1] else x[1]) == 0;
