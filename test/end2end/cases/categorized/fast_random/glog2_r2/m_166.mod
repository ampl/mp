var x{1..3} integer >= -2 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: (x[3] >= -1) or (x[2] != -1) or (x[1] == -2);
