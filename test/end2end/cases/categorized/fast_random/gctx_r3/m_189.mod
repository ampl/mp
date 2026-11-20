var x{1..3} integer >= -1 <= 2;

maximize o: x[1] + 6*x[2] + 5*x[3];

s.t. c1: (if min(x[1], x[2]) != 2 then x[2] else x[2]) >= 0;
s.t. c2: (min(x[1], x[2]) != 2) <==> (min(x[1], x[2]) != 2);
s.t. c3: (max(x[3], x[1]) <= -1) ==> (max(x[3], x[1]) <= -1);
s.t. c4: (x[1] != 0) or (x[1] != 0);
