var x{1..3} integer >= 0 <= 4;
minimize o: (-5)*x[1] + (6)*x[2] + (-1)*x[3];
s.t. c1: (x[1] > 4) ==> (x[2] <= 3) else (x[2] <= 2);
