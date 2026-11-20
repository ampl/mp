var x{1..6} integer >= 0 <= 2;

maximize o: -2*x[1] + 6*x[2] - x[3] - 4*x[4] - x[5] + x[6];

s.t. c1: abs(x[2] - x[1]) == 1;
s.t. c2: (x[6] < 0) ==> (x[5] < 1) else (x[6] != 2);
s.t. c3: ((if x[1] > 2 then 1 else 0) + (if x[6] >= 2 then 1 else 0)) >= 1;
s.t. c4: x[1] == 0;
s.t. c5: x[2] == 1;
