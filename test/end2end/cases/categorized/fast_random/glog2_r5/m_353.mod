var x{1..5} integer >= -2 <= 0;

minimize o: -4*x[1] + 4*x[2] - x[3] - x[4] - 5*x[5];

s.t. c1: ((x[1] == 0) ==> (x[4] > -2)) and ((x[4] > -2) ==> (x[2] >= 0));
s.t. c2: alldiff(x[4], x[1], x[2], x[3]);
s.t. c3: exists {i in 1..5} x[i] == 0;
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
s.t. c5: x[2] == -2;
