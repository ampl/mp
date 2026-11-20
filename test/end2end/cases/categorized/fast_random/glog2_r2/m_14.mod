var x{1..3} integer >= -3 <= 2;

maximize o: 3*x[1] + x[2] + 5*x[3];

s.t. c1: (x[3] < -2) ==> (x[1] == 0) else (x[1] > -1);
