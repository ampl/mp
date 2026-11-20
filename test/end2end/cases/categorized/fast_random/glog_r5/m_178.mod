var x{1..6} integer >= 0 <= 4;

minimize o: 0*x[1] - x[2] - 2*x[3] + 5*x[4] + 4*x[5] - 3*x[6];

s.t. c1: atmost 3 {i in 1..6} (x[i] >= 4);
s.t. c2: not ((x[1] <= 2) and (x[4] <= 4));
s.t. c3: (if x[1] < 0 then x[4] else x[1]) < 3;
s.t. c4: ((x[4] == 4) ==> (x[6] >= 0)) and ((x[6] >= 0) ==> (x[3] >= 4));
