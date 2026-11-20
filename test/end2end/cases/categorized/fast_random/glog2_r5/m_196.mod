var x{1..3} integer >= -2 <= 0;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: (if x[3] <= -2 then x[3] else x[1]) > -1;
s.t. c2: (if x[2] == 0 then x[2] else x[1]) == -2;
s.t. c3: (x[1] > -1) and (x[1] == 0) ==> (x[2] >= -1);
