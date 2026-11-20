var x{1..3} integer >= 0 <= 3;

maximize o: 5*x[1] + 2*x[2] + x[3];

s.t. c1: exists {i in 1..3} x[i] == 1;
s.t. c2: (x[1] == 1) ==> (x[2] >= 1);
s.t. c3: exists {i in 1..3} x[i] == 3;
s.t. c4: ((x[1] > 0) ==> (x[2] >= 2)) and ((x[2] >= 2) ==> (x[2] >= 0));
