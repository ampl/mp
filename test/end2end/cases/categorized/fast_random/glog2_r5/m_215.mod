var x{1..3} integer >= 0 <= 3;

minimize o: 2*x[1] + 0*x[2] + 4*x[3];

s.t. c1: (x[1] >= 1) and (x[2] >= 3) ==> (x[3] >= 1);
s.t. c2: ((if x[3] == 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c3: (if x[1] == 3 then x[3] else x[3]) != 1;
s.t. c4: exists {i in 1..3} x[i] == 3;
s.t. c5: (if x[1] != 3 then x[3] else x[3]) == 1;
