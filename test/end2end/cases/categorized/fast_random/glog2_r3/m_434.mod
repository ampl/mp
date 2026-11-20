var x{1..5} integer >= -2 <= 2;

minimize o: -2*x[1] - 4*x[2] + 4*x[3] + 3*x[4] + 5*x[5];

s.t. c1: (x[4] == -2) <==> (x[1] + x[2] < -2);
s.t. c2: (x[5] != 0) ==> (x[4] < 1) else (x[2] == -1);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: exists {i in 1..5} x[i] == 1;
