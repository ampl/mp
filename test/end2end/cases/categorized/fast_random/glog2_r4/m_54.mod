var x{1..3} integer >= -2 <= 3;

maximize o: 3*x[1] - 5*x[2] + x[3];

s.t. c1: (x[2] >= 3) <==> (x[2] + x[2] != 0);
