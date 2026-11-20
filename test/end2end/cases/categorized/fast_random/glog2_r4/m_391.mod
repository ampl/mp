var x{1..6} integer >= 0 <= 4;

maximize o: 4*x[1] + 4*x[2] + 4*x[3] - 4*x[4] - x[5] + 2*x[6];

s.t. c1: (x[4] == 4) <==> (x[6] + x[4] == 0);
s.t. c2: numberof 4 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 1;
s.t. c3: alldiff(x[3], x[2], x[6]);
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
