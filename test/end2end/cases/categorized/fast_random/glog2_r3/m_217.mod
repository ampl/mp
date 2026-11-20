var x{1..6} integer >= -2 <= 1;
maximize o: (-666)*x[1] + (34)*x[2] + (-2000)*x[3] + (-666)*x[4] + (-666)*x[5] + (-333)*x[6];
s.t. c1: atleast 3 {i in 1..6} (x[i] >= -2);
s.t. c2: (if x[5] > 0 then x[6] else x[5]) != 1;
s.t. c3: atmost 2 {i in 1..6} (x[i] >= 1);
s.t. c4: (x[2] == -1) or (x[5] < -2) or (x[4] <= 1);
