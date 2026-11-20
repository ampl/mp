var x{1..5} integer >= 0 <= 3;

maximize o: 2*x[1] + 2*x[2] - 2*x[3] - 3*x[4] + 2*x[5];

s.t. c1: (x[2] < 2) <==> (x[2] + x[2] < 0);
