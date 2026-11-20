var x{1..3} integer >= 0 <= 2;

maximize o: 4*x[1] - x[2] + 3*x[3];

s.t. c1: numberof 1 in (x[1], x[2], x[3]) >= 1;
s.t. c2: not ((x[3] > 2) and (x[3] > 1));
s.t. c3: abs(x[2] - x[2]) <= 1;
s.t. c4: atmost 1 {i in 1..3} (x[i] >= 2);
s.t. c5: (x[3] >= 1) ==> (x[2] != 2) else (x[3] == 0);
s.t. c6: (x[2] != 0) <==> (x[2] + x[2] > 0);
s.t. c7: alldiff(x[2], x[3]);
