var x{1..5} integer >= 0 <= 2;

maximize o: 5*x[1] - 3*x[2] + 6*x[3] - x[4] + 2*x[5];

s.t. c1: max(x[4], x[1]) - min(x[3], x[5]) >= 3;
s.t. c2: (x[2] < 2) ==> (x[3] != 2);
s.t. c3: alldiff(x[5], x[3]);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
