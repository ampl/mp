var x{1..3} integer >= 0 <= 4;

maximize o: -3*x[1] + 3*x[2] + x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: exists {i in 1..3} x[i] == 3;
s.t. c3: atmost 3 {i in 1..3} (x[i] >= 4);
s.t. c4: not ((x[2] == 3) and (x[2] != 0));
