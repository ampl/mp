var x{1..5} integer >= -2 <= 0;

maximize o: -3*x[1] + x[2] - x[3] + 2*x[4] + 4*x[5];

s.t. c1: ((if x[2] >= -1 then 1 else 0) + (if x[4] >= -1 then 1 else 0)) >= 1;
s.t. c2: abs(x[2] - x[2]) == 2;
s.t. c3: ((x[2] != -2) ==> (x[4] >= -1)) and ((x[4] >= -1) ==> (x[4] >= -2));
s.t. c4: exists {i in 1..5} x[i] == -2;
