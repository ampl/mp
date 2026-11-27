var x{1..4} integer >= -2 <= 0;

maximize o: 6*x[1] - x[2] + 5*x[3] + 6*x[4];

s.t. c1: (x[1] > -3) and (x[3] <= -2) ==> (x[4] >= -3);
