var x{1..3} integer >= 0 <= 2;

maximize o: -5*x[1] + 5*x[2] + 2*x[3] - 6*(if x[1] != x[1] then 1 else 0);

s.t. c1: (x[1] != x[1]) <==> (x[1] == 2);
s.t. c2: not (x[3] < x[3]);
s.t. c3: (x[1] != x[1]) ==> ((abs(x[3] - x[3]) >= 0) or (x[1] == 2));
