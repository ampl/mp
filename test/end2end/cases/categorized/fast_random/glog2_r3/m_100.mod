var x{1..5} integer >= 0 <= 3;

maximize o: 2*x[1] - 3*x[2] + 0*x[3] - 4*x[4] - 4*x[5];

s.t. c1: (if x[3] == 3 then x[1] else x[2]) >= 1;
