var x{1..4} integer >= 0 <= 2;

maximize o: 2000*x[1] + 3000*x[2] + 1000*x[3] - 666*x[4];

s.t. c1: (x[1] < 0) or (x[1] == 2) or (x[1] > 1);
s.t. c2: (x[3] != 0) ==> (x[4] <= 2);
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: abs(x[4] - x[3]) == 1;
