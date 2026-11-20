var x{1..5} integer >= 0 <= 3;

maximize o: 5*x[1] - 2*x[2] - 5*x[3] - x[4] - 5*x[5];

s.t. c1: ((x[5] != 3) ==> (x[3] == 2)) and ((x[3] == 2) ==> (x[5] >= 3));
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) < 3;
s.t. c3: (x[5] >= 3) ==> (x[5] != 0);
s.t. c4: exists {i in 1..5} x[i] == 0;
