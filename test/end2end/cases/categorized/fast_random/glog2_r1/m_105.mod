var x{1..5} integer >= 0 <= 4;

minimize o: 5*x[1] + 2*x[2] + 5*x[3] - 5*x[4] - 3*x[5];

s.t. c1: alldiff(x[3], x[5]);
s.t. c2: (x[4] != 2) ==> (x[4] < 4);
s.t. c3: alldiff(x[1], x[3], x[4]);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
s.t. c5: ((if x[4] <= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
