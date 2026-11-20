var x{1..6} integer >= 0 <= 3;

maximize o: 2*x[1] + 5*x[2] - x[3] - 3*x[4] + 4*x[5] - 5*x[6];

s.t. c1: atleast 1 {i in 1..6} (x[i] >= 1);
s.t. c2: exists {i in 1..6} x[i] == 1;
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 3;
s.t. c4: (x[1] > 3) ==> (x[1] == 0);
