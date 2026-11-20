var x{1..4} integer >= 0 <= 2;

minimize o: -x[1] + 0*x[2] - 5*x[3] - 5*x[4];

s.t. c1: atleast 2 {i in 1..4} (x[i] >= 2);
s.t. c2: (if x[2] < 1 then x[4] else x[2]) >= 2;
s.t. c3: alldiff(x[2], x[3], x[1], x[4]);
s.t. c4: (x[3] < 1) and (x[3] < 2) ==> (x[4] >= 1);
