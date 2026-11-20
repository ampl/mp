var x{1..4} integer >= 0 <= 2;

maximize o: -4*x[1] - 2*x[2] - 5*x[3] - 2*x[4];

s.t. c1: ((if x[1] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: exists {i in 1..4} x[i] == 2;
s.t. c3: (if x[4] > 0 then x[2] else x[1]) >= 2;
s.t. c4: (if x[1] <= 0 then x[1] else x[1]) <= 1;
s.t. c5: ((x[1] == 2) ==> (x[3] == 2)) and ((x[3] == 2) ==> (x[2] >= 2));
s.t. c6: exists {i in 1..4} x[i] == 2;
s.t. c7: forall {i in 1..3} x[i] != x[i+1];
