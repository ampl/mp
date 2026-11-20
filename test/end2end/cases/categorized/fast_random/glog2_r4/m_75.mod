var x{1..6} integer >= 0 <= 4;

minimize o: x[1] + x[2] - 2*x[3] + 3*x[4] + 2*x[5] + x[6];

s.t. c1: abs(x[4] - x[4]) != 0;
s.t. c2: ((if x[3] == 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c3: (x[4] < 3) and (x[4] >= 0) ==> (x[5] >= 3);
s.t. c4: (x[1] < 0) ==> (x[5] <= 2);
s.t. c5: atleast 3 {i in 1..6} (x[i] >= 1);
