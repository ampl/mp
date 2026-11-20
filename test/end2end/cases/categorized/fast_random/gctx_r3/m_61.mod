var x{1..5} integer >= 0 <= 2;

maximize o: 2*x[1] + x[2] + 6*x[3] - 4*x[4] - 3*x[5];

s.t. c1: (x[3] != 1) <==> (x[3] != 1);
s.t. c2: (if x[3] != 1 then x[5] else x[5]) >= 0;
s.t. c3: (x[3] != 1) ==> (x[3] < x[2]);
