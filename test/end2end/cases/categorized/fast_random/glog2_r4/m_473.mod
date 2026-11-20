var x{1..4} integer >= 0 <= 3;

maximize o: -2*x[1] - 3*x[2] - 5*x[3] - 4*x[4];

s.t. c1: atleast 3 {i in 1..4} (x[i] >= 3);
s.t. c2: alldiff(x[2], x[4], x[1]);
s.t. c3: numberof 3 in (x[1], x[2], x[3], x[4]) >= 3;
s.t. c4: x[3] == 3;
s.t. c5: x[1] == 2;
