var x{1..4} integer >= 0 <= 4;

minimize o: -5*x[1] - 5*x[2] + 4*x[3] - 5*x[4];

s.t. c1: alldiff(x[1], x[3]);
s.t. c2: (if x[4] >= 0 then x[4] else x[4]) >= 2;
s.t. c3: alldiff(x[3], x[4], x[1], x[2]);
s.t. c4: (if x[2] >= 1 then x[1] else x[3]) >= 0;
s.t. c5: atleast 1 {i in 1..4} (x[i] >= 1);
