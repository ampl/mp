var x{1..3} integer >= 0 <= 2;

minimize o: 3*x[1] - x[2] - 3*x[3];

s.t. c1: exists {i in 1..3} x[i] == 0;
s.t. c2: (x[2] <= 2) and (x[3] > 2) ==> (x[2] >= 2);
s.t. c3: atmost 2 {i in 1..3} (x[i] >= 2);
s.t. c4: (if x[3] >= 0 then x[1] else x[3]) == 2;
