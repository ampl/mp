var x{1..4} integer >= 0 <= 2;

maximize o: -5*x[1] - 5*x[2] + 4*x[3] - 2*x[4];

s.t. c1: max(x[4], x[4]) - min(x[4], x[4]) == 1;
s.t. c2: atmost 1 {i in 1..4} (x[i] >= 0);
s.t. c3: ((x[3] < 1) ==> (x[1] == 2)) and ((x[1] == 2) ==> (x[3] >= 1));
s.t. c4: exists {i in 1..4} x[i] == 1;
