var x{1..4} integer >= -2 <= 2;

minimize o: 4*x[1] - 5*x[2] + 2*x[3] + 4*x[4];

s.t. c1: alldiff(x[2], x[3]);
s.t. c2: (x[3] > 2) and (x[1] < -2) ==> (x[4] >= 2);
s.t. c3: ((if x[1] == -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c4: alldiff(x[4], x[2]);
s.t. c5: not ((x[4] == 0) and (x[4] != 1));
