var x{1..3} integer >= -2 <= 2;

maximize o: -3*x[1] - 5*x[2] + 5*x[3];

s.t. c1: atmost 2 {i in 1..3} (x[i] >= -2);
s.t. c2: ((x[2] < 0) ==> (x[1] == 2)) and ((x[1] == 2) ==> (x[3] >= 0));
s.t. c3: not ((x[2] > 0) and (x[3] >= -1));
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
