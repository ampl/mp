var x{1..3} integer >= 0 <= 4;

minimize o: 4*x[1] + 6*x[2] + 6*x[3];

s.t. c1: ((x[3] <= 2) ==> (x[3] < 4)) and ((x[3] < 4) ==> (x[2] >= 2));
s.t. c2: alldiff(x[1], x[3], x[2]);
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: x[3] == 1;
s.t. c5: x[1] == 1;
