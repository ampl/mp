var x{1..4} integer >= -2 <= 0;

minimize o: -3*x[1] + 5*x[2] + x[3] - 5*x[4];

s.t. c1: (x[4] <= -2) ==> (x[3] > -2);
s.t. c2: ((x[1] == 0) ==> (x[1] >= 0)) and ((x[1] >= 0) ==> (x[1] >= 0));
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: (x[1] != -1) ==> (x[4] >= -2);
s.t. c5: alldiff(x[2], x[3]);
s.t. c6: numberof 0 in (x[1], x[2], x[3], x[4]) < 3;
