var x{1..5} integer >= 0 <= 4;

maximize o: -x[1] + 4*x[2] + 6*x[3] + 3*x[4] - x[5];

s.t. c1: (if x[2] < 0 then x[4] else x[4]) <= 2;
s.t. c2: (x[3] != 3) <==> (x[2] + x[2] > 2);
s.t. c3: (x[4] > 1) ==> (x[5] == 2);
s.t. c4: atleast 1 {i in 1..5} (x[i] >= 4);
s.t. c5: x[3] == 3;
s.t. c6: x[1] == 4;
