var x{1..4} integer >= 0 <= 4;

maximize o: -2*x[1] - 4*x[2] + 3*x[3] - 5*x[4];

s.t. c1: (x[1] >= 0) ==> (x[1] == 1);
s.t. c2: abs(x[1] - x[1]) == 1;
s.t. c3: (if x[4] >= 2 then x[4] else x[2]) != 0;
s.t. c4: atmost 1 {i in 1..4} (x[i] >= 2);
