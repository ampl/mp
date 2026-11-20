var x{1..3} integer >= -3 <= 1;

minimize o: 3*x[1] - 5*x[2] - 5*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: ((x[2] > -3) ==> (x[3] != 1)) and ((x[3] != 1) ==> (x[1] >= -3));
s.t. c3: ((if x[1] >= -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
