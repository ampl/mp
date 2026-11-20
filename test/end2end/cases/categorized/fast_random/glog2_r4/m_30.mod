var x{1..3} integer >= 0 <= 2;
minimize o: (0)*x[1] + (-3)*x[2] + (-5)*x[3];
s.t. c1: (x[1] != 0) and (x[2] >= 0) ==> (x[3] >= 0);
