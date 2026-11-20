var x{1..5} integer >= -2 <= 1;

minimize o: 6*x[1] + 5*x[2] + 3*x[3] + x[4] + 3*x[5];

s.t. c1: (x[2] != 1) ==> (x[1] <= -1);
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) != 1;
s.t. c3: atmost 2 {i in 1..5} (x[i] >= -1);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
s.t. c5: x[1] == -2;
