var x{1..4} integer >= -2 <= 1;

maximize o: -4*x[1] + 6*x[2] + 0*x[3] + 3*x[4];

s.t. c1: (x[4] <= -1) ==> (x[3] <= 1);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: ((if x[4] >= 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c4: exists {i in 1..4} x[i] == -2;
s.t. c5: atleast 2 {i in 1..4} (x[i] >= 0);
s.t. c6: max(x[1], x[4]) - min(x[1], x[2]) == 1;
s.t. c7: atleast 1 {i in 1..4} (x[i] >= 0);
