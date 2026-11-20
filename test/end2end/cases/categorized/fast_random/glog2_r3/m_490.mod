var x{1..3} integer >= -2 <= 0;

maximize o: 4*x[1] - 3*x[2] + x[3];

s.t. c1: max(x[1], x[2]) - min(x[3], x[1]) < -2;
s.t. c2: (x[2] > 0) ==> (x[2] < 0);
