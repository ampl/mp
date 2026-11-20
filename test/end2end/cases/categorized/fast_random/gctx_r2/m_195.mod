var x{1..5} integer >= 0 <= 2;

maximize o: 0*x[1] + 4*x[2] + 2*x[3] + 0*x[4] + x[5];

s.t. c1: (min(x[2], x[3]) == 1) <==> (x[1] >= 2);
s.t. c2: not (min(x[2], x[3]) == 1);
s.t. c3: (if min(x[2], x[3]) == 1 then x[4] else x[4]) < 2;
s.t. c4: (if min(x[2], x[3]) == 1 then x[5] else x[5]) <= 1;
