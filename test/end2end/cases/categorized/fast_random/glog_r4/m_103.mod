var x{1..4} integer >= 0 <= 2;

maximize o: 2*x[1] - 5*x[2] + 3*x[3] - 3*x[4];

s.t. c1: atleast 2 {i in 1..4} (x[i] >= 1);
s.t. c2: exists {i in 1..4} x[i] == 1;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
