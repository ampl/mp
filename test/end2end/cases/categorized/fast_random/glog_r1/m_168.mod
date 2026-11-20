var x{1..3} integer >= 0 <= 4;

minimize o: 6*x[1] - 5*x[2] - 3*x[3];

s.t. c1: ((x[1] >= 2) ==> (x[2] > 3)) and ((x[2] > 3) ==> (x[2] >= 2));
s.t. c2: (if x[2] >= 4 then x[2] else x[1]) <= 2;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: ((if x[2] != 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
