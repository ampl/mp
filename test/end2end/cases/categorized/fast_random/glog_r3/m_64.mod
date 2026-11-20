var x{1..5} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (-3)*x[2] + (-2)*x[3] + (-2)*x[4] + (1)*x[5];
s.t. c1: (x[1] >= 0) ==> (x[1] > 0);
s.t. c2: exists {i in 1..5} x[i] == 0;
s.t. c3: (if x[5] != 2 then x[2] else x[3]) != 1;
s.t. c4: atmost 1 {i in 1..5} (x[i] >= 2);
