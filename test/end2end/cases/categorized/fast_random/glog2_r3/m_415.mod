var x{1..4} integer >= 0 <= 3;

maximize o: 4*x[1] - 3*x[2] + 6*x[3] - x[4];

s.t. c1: max(x[2], x[2]) - min(x[2], x[2]) <= 2;
s.t. c2: ((x[2] > 3) ==> (x[2] == 3)) and ((x[2] == 3) ==> (x[1] >= 3));
s.t. c3: atleast 2 {i in 1..4} (x[i] >= 0);
s.t. c4: (if x[3] >= 2 then x[1] else x[2]) <= 2;
s.t. c5: numberof 0 in (x[1], x[2], x[3], x[4]) < 1;
s.t. c6: (if x[3] != 2 then x[2] else x[4]) <= 1;
s.t. c7: ((if x[2] != 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
