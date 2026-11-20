var x{1..3} integer >= 0 <= 4;

maximize o: -x[1] + 6*x[2] - 5*x[3];

s.t. c1: exists {i in 1..3} x[i] == 3;
s.t. c2: ((x[2] < 2) ==> (x[2] > 0)) and ((x[2] > 0) ==> (x[3] >= 2));
s.t. c3: ((if x[1] > 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c4: abs(x[3] - x[3]) != 2;
