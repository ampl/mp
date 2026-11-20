var x{1..5} integer >= -2 <= 2;

maximize o: 2*x[1] + x[2] + 4*x[3] + 5*x[4] - 3*x[5];

s.t. c1: (if x[2] != 0 then x[3] else x[3]) > 2;
s.t. c2: (x[5] < 2) and (x[4] <= 2) ==> (x[4] >= 2);
s.t. c3: (x[4] > -2) <==> (x[5] + x[5] < -2);
s.t. c4: x[1] == 1;
