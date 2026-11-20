var x{1..6} integer >= -2 <= 1;

maximize o: 5*x[1] + 3*x[2] + 3*x[3] - 3*x[4] + 0*x[5] - 4*x[6];

s.t. c1: (x[3] > 1) and (x[5] >= 0) ==> (x[6] >= 1);
s.t. c2: not ((x[2] < 0) and (x[3] > 1));
s.t. c3: (x[1] <= -1) ==> (x[3] <= -2) else (x[5] > -2);
s.t. c4: atleast 2 {i in 1..6} (x[i] >= -1);
s.t. c5: x[1] == -2;
