var x{1..4} integer >= 0 <= 2;

minimize o: -666*x[1] - 3000*x[2] + 3000*x[3] + 34*x[4];

s.t. c1: (if x[3] <= 0 then x[2] else x[3]) <= 1;
s.t. c2: atmost 1 {i in 1..4} (x[i] >= 0);
