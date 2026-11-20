var x{1..3} integer >= 0 <= 2;

maximize o: -5*x[1] + x[2] + 2*x[3];

s.t. c1: numberof 2 in (x[1], x[2], x[3]) >= 0;
s.t. c2: atmost 2 {i in 1..3} (x[i] >= 0);
s.t. c3: (x[1] != 1) ==> (x[2] > 1) else (x[3] != 0);
