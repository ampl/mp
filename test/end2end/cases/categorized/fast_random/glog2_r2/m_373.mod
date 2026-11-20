var x{1..3} integer >= 0 <= 3;

maximize o: 6*x[1] - 4*x[2] + 0*x[3];

s.t. c1: (x[2] >= 2) <==> (x[3] + x[3] > 1);
s.t. c2: alldiff(x[2], x[1], x[3]);
s.t. c3: atleast 2 {i in 1..3} (x[i] >= 1);
s.t. c4: x[2] == 2;
s.t. c5: x[1] == 0;
