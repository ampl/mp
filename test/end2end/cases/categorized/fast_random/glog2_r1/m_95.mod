var x{1..6} integer >= -2 <= 0;

minimize o: -5*x[1] - 5*x[2] - 3*x[3] + 3*x[4] - 3*x[5] + 3*x[6];

s.t. c1: ((x[4] > -2) ==> (x[6] > 0)) and ((x[6] > 0) ==> (x[5] >= -2));
s.t. c2: ((x[2] >= -1) ==> (x[3] == -1)) and ((x[3] == -1) ==> (x[4] >= -1));
s.t. c3: (if x[1] != 0 then x[6] else x[5]) >= -1;
s.t. c4: (if x[6] != 0 then x[5] else x[1]) > -2;
s.t. c5: max(x[5], x[5]) - min(x[2], x[4]) < -2;
s.t. c6: ((if x[4] < 0 then 1 else 0) + (if x[6] >= 0 then 1 else 0)) >= 1;
s.t. c7: atleast 1 {i in 1..6} (x[i] >= -1);
