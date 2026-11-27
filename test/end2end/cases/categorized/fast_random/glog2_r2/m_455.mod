var x{1..3} integer >= 0 <= 2;

maximize o: 2*x[1] + 6*x[2] - 2*x[3];

s.t. c1: alldiff(x[1], x[2], x[3]);
s.t. c2: ((x[3] != 2) ==> (x[1] < 0)) and ((x[1] < 0) ==> (x[1] >= 2));
s.t. c3: abs(x[1] - x[1]) != 0;
s.t. c4: alldiff(x[2], x[3], x[1]);
s.t. c5: not ((x[2] < 2) and (x[1] > 2));
