var x{1..6} integer >= 0 <= 3;

maximize o: -x[1] + 6*x[2] + 3*x[3] + 2*x[4] + 5*x[5] + 4*x[6];

s.t. c1: exists {i in 1..6} x[i] == 0;
s.t. c2: (if x[6] > 2 then x[1] else x[4]) < 2;
s.t. c3: atmost 3 {i in 1..6} (x[i] >= 3);
s.t. c4: alldiff(x[5], x[2], x[4], x[1]);
s.t. c5: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 2;
s.t. c6: (x[5] == 0) and (x[2] < 1) ==> (x[4] >= 0);
s.t. c7: abs(x[2] - x[4]) != 1;
