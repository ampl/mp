var x{1..4} integer >= 0 <= 4;

maximize o: -666*x[1] + 51*x[2] - 333*x[3] + 1000*x[4];

s.t. c1: (x[4] != 2) ==> (x[1] >= 1) else (x[4] > 1);
