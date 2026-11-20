var x{1..3} integer >= -2 <= 2;

maximize o: 2000*x[1] - 2000*x[2] + 2000*x[3];

s.t. c1: exists {i in 1..3} x[i] == -1;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: exists {i in 1..3} x[i] == 0;
s.t. c4: atleast 3 {i in 1..3} (x[i] >= -2);
