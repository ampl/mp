var x{1..6} integer >= 0 <= 2;

minimize o: 6*x[1] - 3*x[2] + 6*x[3] - 3*x[4] + 4*x[5] + 4*x[6];

s.t. c1: ((x[6] != 0) ==> (x[4] <= 0)) and ((x[4] <= 0) ==> (x[3] >= 0));
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 3;
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
s.t. c4: ((x[3] == 1) ==> (x[6] > 0)) and ((x[6] > 0) ==> (x[5] >= 1));
s.t. c5: forall {i in 1..5} x[i] != x[i+1];
s.t. c6: max(x[6], x[3]) - min(x[5], x[3]) < 2;
s.t. c7: alldiff(x[1], x[4]);
