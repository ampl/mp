var x{1..3} integer >= 0 <= 2;

minimize o: -2*x[1] - x[2] + 6*x[3];

s.t. c1: not ((x[3] < 1) and (x[3] == 1));
s.t. c2: atleast 2 {i in 1..3} (x[i] >= 0);
s.t. c3: ((if x[1] == 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c4: (x[1] != 1) ==> (x[1] <= 0) else (x[2] == 0);
