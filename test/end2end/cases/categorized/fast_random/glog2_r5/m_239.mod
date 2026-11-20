var x{1..6} integer >= -2 <= 0;

maximize o: -5*x[1] + 4*x[2] + 6*x[3] - 4*x[4] - 5*x[5] + 6*x[6];

s.t. c1: exists {i in 1..6} x[i] == 0;
s.t. c2: ((x[5] >= -1) ==> (x[6] == 0)) and ((x[6] == 0) ==> (x[1] >= -1));
s.t. c3: ((if x[2] < 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c4: (if x[1] < -1 then x[1] else x[1]) == -2;
s.t. c5: exists {i in 1..6} x[i] == 0;
s.t. c6: (if x[1] < -1 then x[1] else x[1]) == -2;
