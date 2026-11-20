var x{1..5} integer >= -2 <= 1;

maximize o: -4*x[1] - 4*x[2] - 2*x[3] - 4*x[4] + 6*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: ((x[2] != -1) ==> (x[4] > 0)) and ((x[4] > 0) ==> (x[1] >= -1));
s.t. c3: atmost 1 {i in 1..5} (x[i] >= -1);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
s.t. c5: ((x[5] < -2) ==> (x[4] >= -2)) and ((x[4] >= -2) ==> (x[1] >= -2));
s.t. c6: (x[1] <= -1) ==> (x[4] == -2) else (x[2] != 1);
s.t. c7: alldiff(x[5], x[4], x[1]);
