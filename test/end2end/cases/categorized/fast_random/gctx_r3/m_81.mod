var x{1..3} integer >= 0 <= 2;

maximize o: -3*x[1] - 4*x[2] - 2*x[3];

s.t. c1: (abs(x[1] - x[1]) == 1) or (max(x[3], x[1]) == 1);
s.t. c2: (abs(x[1] - x[1]) == 1) and ((x[3] >= 1) or not (x[3] >= 1));
