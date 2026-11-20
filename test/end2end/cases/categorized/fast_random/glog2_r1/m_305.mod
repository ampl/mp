var x{1..3} integer >= 0 <= 4;

minimize o: 5*x[1] + 3*x[2] + x[3];

s.t. c1: (x[1] != 1) ==> (x[3] <= 4) else (x[3] > 4);
s.t. c2: alldiff(x[1], x[3]);
s.t. c3: (x[2] > 3) and (x[2] != 0) ==> (x[3] >= 3);
s.t. c4: alldiff(x[2], x[1]);
s.t. c5: atmost 1 {i in 1..3} (x[i] >= 1);
s.t. c6: (x[1] > 0) <==> (x[3] + x[2] > 4);
s.t. c7: max(x[1], x[1]) - min(x[3], x[3]) > 4;
