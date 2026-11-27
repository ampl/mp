var x{1..6} integer >= -2 <= 1;

maximize o: 3*x[1] + 6*x[2] - 2*x[3] + 4*x[4] + 5*x[5] + 5*x[6];

s.t. c1: not ((x[1] < -2) and (x[2] > -2));
s.t. c2: numberof -3 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 3;
s.t. c3: atmost 3 {i in 1..6} (x[i] >= 1);
s.t. c4: alldiff(x[2], x[3]);
