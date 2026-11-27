var x{1..4} integer >= -2 <= 1;

maximize o: -x[1] - x[2] - 4*x[3] + 3*x[4];

s.t. c1: (x[3] >= -1) <==> (x[1] + x[2] == 0);
