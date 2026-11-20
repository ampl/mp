var x{1..3} integer >= 0 <= 2;

maximize o: -x[1] - 3*x[2] - 3*x[3];

s.t. c1: (x[2] > -1) <==> (x[2] + x[2] <= 2);
