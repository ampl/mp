var x{1..3} integer >= 0 <= 3;

maximize o: -2*x[1] - x[2] - 3*x[3];

s.t. c1: (x[3] != 3) ==> (max(x[3], x[2]) > 2);
s.t. c2: (x[3] != 3) ==> (abs(x[1] - x[2]) < 1) else (max(x[3], x[2]) > 2);
