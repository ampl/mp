var x{1..4} integer >= -2 <= 2;

maximize o: -2*x[1] + x[2] + 5*x[3] + 5*x[4];

s.t. c1: max(x[3], x[3]) - min(x[4], x[4]) < 0;
s.t. c2: not ((x[1] >= 0) and (x[3] > 0));
s.t. c3: alldiff(x[3], x[2]);
s.t. c4: (if x[3] > -2 then x[3] else x[4]) >= 1;
s.t. c5: (x[4] < 2) or (x[2] > -1) or (x[2] <= -1);
s.t. c6: (x[2] == 0) <==> (x[4] + x[4] != 1);
s.t. c7: forall {i in 1..3} x[i] != x[i+1];
