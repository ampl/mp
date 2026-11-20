var x{1..5} integer >= -2 <= 2;

maximize o: -4*x[1] + x[2] + 2*x[3] - 4*x[4] + 5*x[5];

s.t. c1: (x[4] >= -3) <==> (x[2] + x[3] < 2);
s.t. c2: ((x[2] > 3) ==> (x[2] != -3)) and ((x[2] != -3) ==> (x[4] >= 3));
s.t. c3: alldiff(x[5], x[2], x[1]);
s.t. c4: atleast 1 {i in 1..5} (x[i] >= 2);
