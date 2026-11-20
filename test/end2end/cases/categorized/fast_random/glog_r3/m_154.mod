var x{1..5} integer >= 0 <= 3;

minimize o: -3*x[1] - 4*x[2] - 4*x[3] + 6*x[4] - x[5];

s.t. c1: atleast 1 {i in 1..5} (x[i] >= 2);
s.t. c2: (if x[5] > 3 then x[3] else x[4]) == 3;
s.t. c3: ((x[5] != 2) ==> (x[4] == 3)) and ((x[4] == 3) ==> (x[3] >= 2));
s.t. c4: ((if x[4] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
