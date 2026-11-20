var x{1..6} integer >= 0 <= 2;

maximize o: -x[1] + 3*x[2] + x[3] + 6*x[4] - x[5] - 3*x[6];

s.t. c1: (if x[6] > 0 then x[2] else x[2]) >= 0;
s.t. c2: (x[1] <= 2) ==> (x[5] <= 2);
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 3;
s.t. c4: (x[2] >= 0) <==> (x[3] + x[3] <= 2);
s.t. c5: atleast 1 {i in 1..6} (x[i] >= 2);
s.t. c6: abs(x[3] - x[3]) <= 2;
s.t. c7: abs(x[3] - x[5]) < 1;
