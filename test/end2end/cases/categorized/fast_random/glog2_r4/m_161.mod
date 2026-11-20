var x{1..4} integer >= 0 <= 3;

minimize o: -x[1] - 5*x[2] + 6*x[3] + 3*x[4];

s.t. c1: ((if x[3] <= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c2: ((x[2] < -1) ==> (x[4] < 3)) and ((x[4] < 3) ==> (x[4] >= -1));
s.t. c3: ((x[1] < 3) ==> (x[1] > 3)) and ((x[1] > 3) ==> (x[4] >= 3));
s.t. c4: exists {i in 1..4} x[i] == 3;
