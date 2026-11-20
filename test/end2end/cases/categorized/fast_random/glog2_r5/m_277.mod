var x{1..4} integer >= 0 <= 3;
minimize o: (-3000)*x[1] + (-999)*x[2] + (-666)*x[3] + (3000)*x[4];
s.t. c1: atmost 3 {i in 1..4} (x[i] >= 0);
s.t. c2: (x[1] <= 0) ==> (x[4] != 3) else (x[2] != 2);
