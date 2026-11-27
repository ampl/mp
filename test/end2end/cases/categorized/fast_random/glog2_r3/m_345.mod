var x{1..5} integer >= -2 <= 2;

minimize o: 2*x[1] + 2*x[2] + 2*x[3] - x[4] + 0*x[5];

s.t. c1: atmost 1 {i in 1..5} (x[i] >= 0);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: abs(x[1] - x[5]) > 0;
s.t. c4: (x[3] <= 2) ==> (x[1] >= -1);
s.t. c5: abs(x[1] - x[2]) > 1;
