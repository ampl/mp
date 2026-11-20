var x{1..4} integer >= 0 <= 2;
minimize o: (3)*x[1] + (3)*x[2] + (1)*x[3] + (-5)*x[4];
s.t. c1: (x[4] != x[4]) ==> (x[4] != x[4]);
s.t. c2: (x[3] > 1) or (abs(x[1] - x[2]) > 0);
s.t. c3: (x[3] > 1) or (x[3] > 1);
s.t. c4: (x[4] != x[4]) and ((abs(x[1] - x[2]) > 0) or not (x[3] > 1));
