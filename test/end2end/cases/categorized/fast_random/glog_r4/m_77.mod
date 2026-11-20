var x{1..4} integer >= 0 <= 3;

minimize o: 4*x[1] + x[2] + 6*x[3] + 4*x[4];

s.t. c1: ((if x[2] == 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[3] < 0) ==> (x[1] == 0);
s.t. c3: atleast 1 {i in 1..4} (x[i] >= 0);
s.t. c4: ((if x[4] > 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
