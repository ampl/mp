var x{1..3} integer >= -2 <= 3;

maximize o: x[1] + 6*x[2] + 5*x[3];

s.t. c1: ((if x[3] >= 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: ((x[1] > 2) ==> (x[2] < 1)) and ((x[2] < 1) ==> (x[2] >= 2));
s.t. c3: abs(x[3] - x[1]) <= 1;
s.t. c4: exists {i in 1..3} x[i] == 0;
