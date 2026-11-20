var x{1..5} integer >= 0 <= 2;

minimize o: -x[1] - 3*x[2] - x[3] + 4*x[4] - 5*x[5];

s.t. c1: ((if x[1] >= 3 then 1 else 0) + (if x[5] >= 3 then 1 else 0)) >= 1;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: atleast 3 {i in 1..5} (x[i] >= 2);
