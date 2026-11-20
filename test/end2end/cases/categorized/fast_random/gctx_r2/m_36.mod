var x{1..4} integer >= -1 <= 1;
minimize o: (4)*x[1] + (4)*x[2] + (-5)*x[3] + (4)*x[4];
s.t. c1: (if x[4] + x[1] != 0 then x[4] else 0) + (if x[2] + x[3] < 1 then x[1] else 0) > 0;
s.t. c2: (x[2] + x[3] < 1) ==> ((x[2] + x[3] < 1) or (x[4] + x[1] != 0));
