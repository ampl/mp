var x{1..5} integer >= 0 <= 4;

minimize o: 2*x[1] - 4*x[2] - 5*x[3] - 4*x[4] + 6*x[5];

s.t. c1: atleast 2 {i in 1..5} (x[i] >= 0);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (x[1] == 0) or (x[3] >= 4) or (x[3] == -1);
s.t. c4: ((x[3] == -1) ==> (x[4] < 0)) and ((x[4] < 0) ==> (x[1] >= -1));
