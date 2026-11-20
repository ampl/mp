var x{1..5} integer >= -2 <= 1;

minimize o: -4*x[1] - 4*x[2] + 5*x[3] + 5*x[4] + 4*x[5];

s.t. c1: (if x[3] >= -2 then x[4] else x[1]) == 1;
s.t. c2: atmost 2 {i in 1..5} (x[i] >= 2);
s.t. c3: exists {i in 1..5} x[i] == 1;
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
