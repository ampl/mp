var x{1..4} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (x[3] == 0) <==> (x[4] + x[4] >= 2);
s.t. c2: abs(x[2] - x[2]) > 2;
s.t. c3: (x[3] != 1) and (x[1] == 2) ==> (x[4] >= 1);
s.t. c4: atleast 3 {i in 1..4} (x[i] >= 0);
