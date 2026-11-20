var x{1..5} integer >= 0 <= 2;

minimize o: x[1] - 4*x[2] + 3*x[3] + 6*x[4] + 4*x[5];

s.t. c1: atleast 3 {i in 1..5} (x[i] >= 0);
s.t. c2: exists {i in 1..5} x[i] == 0;
s.t. c3: abs(x[3] - x[3]) < 2;
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
s.t. c5: alldiff(x[4], x[5], x[1]);
