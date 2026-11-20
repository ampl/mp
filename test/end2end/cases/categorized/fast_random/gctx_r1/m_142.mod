var x{1..4} integer >= 0 <= 3;

maximize o: 4*x[1] - x[2] + 3*x[3] + 3*x[4];

s.t. c1: not (x[4] < 3);
s.t. c2: (abs(x[2] - x[2]) <= 1) ==> (x[4] < 3) else (abs(x[2] - x[2]) <= 1);
