var x{1..3} integer >= 0 <= 2;

maximize o: 2*x[1] + 5*x[2] + 4*x[3] - (if x[1] != 1 then 1 else 0);

s.t. c1: not ((x[1] < 0) and (x[1] < 0));
s.t. c2: (x[1] < 0) ==> ((x[1] != 1) or (x[1] != 1));
s.t. c3: (if x[1] <= 1 then x[2] else 0) + (if x[1] <= 1 then x[1] else 0) < 0;
s.t. c4: (x[3] == x[1]) <==> (x[1] < 0);
