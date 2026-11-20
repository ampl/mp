var x{1..6} integer >= 0 <= 2;

maximize o: 5*x[1] + 3*x[2] - 5*x[3] + 2*x[4] + 5*x[5] + 3*x[6];

s.t. c1: (x[5] == 0) or (x[4] == 0) or (x[2] > 2);
s.t. c2: alldiff(x[1], x[4], x[6]);
s.t. c3: exists {i in 1..6} x[i] == 1;
s.t. c4: (if x[3] < 2 then x[6] else x[3]) > 1;
s.t. c5: forall {i in 1..5} x[i] != x[i+1];
s.t. c6: (x[3] != 0) or (x[6] > 0) or (x[3] > 2);
