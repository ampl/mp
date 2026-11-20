var x{1..3} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (-2)*x[2] + (3)*x[3];
s.t. c1: not ((abs(x[1] - x[2]) <= 0) and (abs(x[1] - x[2]) <= 0));
s.t. c2: (abs(x[1] - x[2]) <= 0) ==> (abs(x[1] - x[2]) <= 0) else (x[2] >= x[2]);
