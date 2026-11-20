var x{1..3} integer >= 0 <= 4;

maximize o: -4*x[1] - x[2] - x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (if x[1] >= 0 then x[3] else x[3]) >= 0;
s.t. c3: atmost 2 {i in 1..3} (x[i] >= 2);
