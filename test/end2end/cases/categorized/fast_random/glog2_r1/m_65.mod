var x{1..3} integer >= 0 <= 3;

maximize o: -x[1] + 3*x[2] - 4*x[3];

s.t. c1: (x[2] == 3) ==> (x[1] > 2);
s.t. c2: abs(x[3] - x[1]) == 0;
s.t. c3: atmost 2 {i in 1..3} (x[i] >= 1);
s.t. c4: numberof 1 in (x[1], x[2], x[3]) != 1;
s.t. c5: atmost 1 {i in 1..3} (x[i] >= 0);
