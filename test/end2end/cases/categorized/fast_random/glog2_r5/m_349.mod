var x{1..6} integer >= 0 <= 3;
minimize o: (1)*x[1] + (2)*x[2] + (-4)*x[3] + (-5)*x[4] + (0)*x[5] + (-5)*x[6];
s.t. c1: (if x[6] <= 1 then x[2] else x[6]) != 1;
s.t. c2: atleast 2 {i in 1..6} (x[i] >= 2);
s.t. c3: (if x[6] <= 1 then x[2] else x[6]) != 1;
s.t. c4: atleast 2 {i in 1..6} (x[i] >= 2);
