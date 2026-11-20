var x{1..5} integer >= 0 <= 2;

minimize o: x[1] + 4*x[2] - 5*x[3] - 3*x[4] - 3*x[5];

s.t. c1: (x[2] == 1) ==> (x[2] < 1);
s.t. c2: abs(x[4] - x[2]) >= 1;
s.t. c3: (if x[1] < 1 then x[3] else x[4]) >= 0;
s.t. c4: not ((x[5] > 1) and (x[3] < 2));
s.t. c5: (x[5] > 2) ==> (x[2] > 1);
s.t. c6: alldiff(x[2], x[4]);
s.t. c7: alldiff(x[1], x[2]);
