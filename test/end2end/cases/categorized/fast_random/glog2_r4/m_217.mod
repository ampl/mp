var x{1..3} integer >= 0 <= 3;

maximize o: 51*x[1] - 333*x[2] - 1000*x[3];

s.t. c1: (x[1] < 2) or (x[1] <= 0) or (x[2] < 0);
