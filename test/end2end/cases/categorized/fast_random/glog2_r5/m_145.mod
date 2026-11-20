var x{1..3} integer >= -2 <= 2;

maximize o: -3*x[1] - 3*x[2] + 6*x[3];

s.t. c1: (x[1] <= 0) ==> (x[2] < 1);
s.t. c2: atleast 2 {i in 1..3} (x[i] >= -2);
s.t. c3: alldiff(x[2], x[3]);
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
s.t. c5: not ((x[1] < 1) and (x[1] >= -1));
