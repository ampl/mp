var x{1..5} integer >= -2 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: atleast 1 {i in 1..5} (x[i] >= -1);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: not ((x[4] != -1) and (x[4] != 2));
s.t. c4: (x[2] > -1) <==> (x[1] + x[5] < 1);
