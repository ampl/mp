var x{1..7} integer >= 0 <= 1;

maximize o: 5*x[1] + x[2] - 3*x[3] + 2*x[4] - 2*x[5] - x[6] + 0*x[7];

s.t. c1: (x[7] > 1) and (x[3] == 0) ==> (x[6] >= 1);
s.t. c2: (x[5] <= 0) ==> (x[1] > 1);
s.t. c3: atleast 1 {i in 1..7} (x[i] >= 1);
s.t. c4: alldiff(x[7], x[5], x[3]);
