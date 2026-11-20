var x{1..5} integer >= -2 <= 0;

maximize o: -5*x[1] - x[2] + x[3] - 2*x[4] + 6*x[5];

s.t. c1: exists {i in 1..5} x[i] == -2;
s.t. c2: (if x[2] != -2 then x[3] else x[1]) >= -2;
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: exists {i in 1..5} x[i] == -2;
s.t. c5: forall {i in 1..4} x[i] != x[i+1];
