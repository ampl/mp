var x{1..3} integer >= -2 <= 2;

minimize o: -x[1] - x[2] - 2*x[3];

s.t. c1: ((x[2] != 1) ==> (x[3] != 2)) and ((x[3] != 2) ==> (x[3] >= 1));
s.t. c2: exists {i in 1..3} x[i] == 0;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: max(x[2], x[1]) - min(x[2], x[2]) >= -1;
s.t. c5: ((x[2] != 1) ==> (x[3] != 2)) and ((x[3] != 2) ==> (x[3] >= 1));
s.t. c6: max(x[2], x[1]) - min(x[2], x[2]) >= -1;
