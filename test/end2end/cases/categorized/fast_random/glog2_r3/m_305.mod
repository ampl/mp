var x{1..5} integer >= 0 <= 3;
minimize o: (1)*x[1] + (4)*x[2] + (2)*x[3] + (2)*x[4] + (-5)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: atmost 3 {i in 1..5} (x[i] >= 3);
s.t. c3: (x[2] > 1) or (x[2] <= 1) or (x[2] < 2);
s.t. c4: atmost 2 {i in 1..5} (x[i] >= 1);
s.t. c5: exists {i in 1..5} x[i] == 2;
s.t. c6: (x[2] <= 3) ==> (x[3] > 0);
