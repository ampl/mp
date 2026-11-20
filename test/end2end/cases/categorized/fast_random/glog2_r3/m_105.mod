var x{1..6} integer >= -2 <= 2;

maximize o: 0*x[1] + 2*x[2] + x[3] - 2*x[4] - 3*x[5] + 6*x[6];

s.t. c1: alldiff(x[4], x[5], x[2]);
s.t. c2: (x[6] > 1) ==> (x[6] == -2);
s.t. c3: max(x[1], x[5]) - min(x[3], x[1]) == 1;
s.t. c4: max(x[5], x[2]) - min(x[1], x[4]) >= -1;
s.t. c5: (if x[6] != -2 then x[3] else x[5]) >= 0;
s.t. c6: abs(x[4] - x[4]) > 3;
s.t. c7: atleast 2 {i in 1..6} (x[i] >= 1);
