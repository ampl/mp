var x{1..6} integer >= -2 <= 1;

minimize o: -x[1] + 6*x[2] - 5*x[3] + 4*x[4] + 4*x[5] + 0*x[6];

s.t. c1: atleast 3 {i in 1..6} (x[i] >= -1);
s.t. c2: (x[1] > -1) and (x[3] <= 0) ==> (x[3] >= -1);
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 3;
s.t. c4: alldiff(x[2], x[4], x[5]);
