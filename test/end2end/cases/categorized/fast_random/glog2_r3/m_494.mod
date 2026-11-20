var x{1..4} integer >= -2 <= 2;

maximize o: -2*x[1] + 2*x[2] - x[3] - 3*x[4];

s.t. c1: atmost 2 {i in 1..4} (x[i] >= -2);
s.t. c2: ((x[3] <= 0) ==> (x[1] != -2)) and ((x[1] != -2) ==> (x[3] >= 0));
