var x{1..3} integer >= 0 <= 2;
minimize o: (-2)*x[1] + (1)*x[2] + (4)*x[3];
s.t. c1: not (min(x[3], x[2]) <= 0);
s.t. c2: (abs(x[3] - x[1]) == 2) ==> ((abs(x[3] - x[1]) == 2) or (abs(x[3] - x[1]) == 2));
