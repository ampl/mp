var x{1..3} integer >= 0 <= 4;

maximize o: -4*x[1] + 4*x[2] - x[3];

s.t. c1: (if x[2] >= 4 then x[1] else x[3]) == 1;
s.t. c2: (x[3] < 3) or (x[2] == 1) or (x[1] < 1);
s.t. c3: alldiff(x[2], x[3], x[1]);
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
