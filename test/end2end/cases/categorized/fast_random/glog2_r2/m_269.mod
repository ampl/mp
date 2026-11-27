var x{1..3} integer >= -2 <= 0;

maximize o: 2*x[1] + 3*x[2] + 5*x[3];

s.t. c1: exists {i in 1..3} x[i] == -2;
s.t. c2: (x[2] < 0) and (x[2] > 0) ==> (x[1] >= 0);
s.t. c3: alldiff(x[1], x[2], x[3]);
s.t. c4: numberof 0 in (x[1], x[2], x[3]) > 3;
s.t. c5: exists {i in 1..3} x[i] == -2;
s.t. c6: numberof 0 in (x[1], x[2], x[3]) > 3;
