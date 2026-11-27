var x{1..4} integer >= 0 <= 3;

minimize o: -2*x[1] + 3*x[2] - 4*x[3] + 4*x[4];

s.t. c1: (x[2] >= 0) ==> (x[1] > 0) else (x[2] > 2);
s.t. c2: atleast 2 {i in 1..4} (x[i] >= 1);
s.t. c3: (if x[3] > 0 then x[2] else x[3]) > 1;
s.t. c4: abs(x[2] - x[1]) >= 1;
s.t. c5: (x[1] < 2) and (x[3] > 1) ==> (x[1] >= 2);
