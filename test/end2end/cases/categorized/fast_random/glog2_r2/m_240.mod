var x{1..5} integer >= 0 <= 4;
maximize o: (6)*x[1] + (-5)*x[2] + (1)*x[3] + (6)*x[4] + (4)*x[5];
s.t. c1: abs(x[4] - x[5]) >= 1;
