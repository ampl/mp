var x{1..4} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (1)*x[2] + (5)*x[3] + (2)*x[4];
s.t. c1: (x[1] >= 2) ==> (x[2] == 1);
s.t. c2: (x[1] != 1) ==> (x[1] != 2);
s.t. c3: (x[1] != 3) and (x[1] == 2) ==> (x[2] >= 3);
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
s.t. c5: (x[4] == 1) ==> (x[1] != 2) else (x[3] < 4);
s.t. c6: atleast 1 {i in 1..4} (x[i] >= 2);
s.t. c7: forall {i in 1..3} x[i] != x[i+1];
