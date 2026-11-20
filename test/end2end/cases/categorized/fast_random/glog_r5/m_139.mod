var x{1..5} integer >= 0 <= 4;

minimize o: -3*x[1] - 2*x[2] - x[3] + 3*x[4] + 3*x[5];

s.t. c1: (if x[4] < 4 then x[3] else x[3]) >= 0;
s.t. c2: (x[4] != 1) and (x[2] >= 2) ==> (x[1] >= 1);
s.t. c3: atmost 2 {i in 1..5} (x[i] >= 4);
