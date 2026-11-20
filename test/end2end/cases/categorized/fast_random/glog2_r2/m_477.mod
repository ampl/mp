var x{1..3} integer >= 0 <= 2;

maximize o: 1000*x[1] + 3000*x[2] + 3000*x[3];

s.t. c1: atmost 1 {i in 1..3} (x[i] >= 0);
s.t. c2: (x[3] > 2) ==> (x[2] <= 1);
s.t. c3: exists {i in 1..3} x[i] == 1;
