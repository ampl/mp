var x{1..4} integer >= 0 <= 3;

minimize o: -x[1] + 6*x[2] + 2*x[3] + 0*x[4];

s.t. c1: exists {i in 1..4} x[i] == 2;
s.t. c2: ((x[4] == 2) ==> (x[3] < 2)) and ((x[3] < 2) ==> (x[4] >= 2));
s.t. c3: ((if x[1] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c4: exists {i in 1..4} x[i] == 1;
