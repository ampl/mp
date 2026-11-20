var x{1..4} integer >= -2 <= 1;

maximize o: 3000*x[1] + 1000*x[2] + 2000*x[3] - 666*x[4];

s.t. c1: not ((x[2] != -2) and (x[4] != 0));
s.t. c2: (x[3] >= 1) <==> (x[4] + x[2] == -1);
s.t. c3: alldiff(x[4], x[1], x[2], x[3]);
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
