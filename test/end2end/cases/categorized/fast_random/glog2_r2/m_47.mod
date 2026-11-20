var x{1..3} integer >= 0 <= 2;
minimize o: (-666)*x[1] + (-333)*x[2] + (51)*x[3];
s.t. c1: atmost 1 {i in 1..3} (x[i] >= 2);
s.t. c2: (x[2] <= 1) and (x[2] < 2) ==> (x[1] >= 1);
s.t. c3: (x[1] >= 0) or (x[2] >= 1) or (x[1] != 1);
