var x{1..3} integer >= -1 <= 2;

maximize o: 6*x[1] + 2*x[2] + 3*x[3] - 5*(if x[3] == 0 then 1 else 0);

s.t. c1: not ((x[3] == 0) and (x[3] == 0));
s.t. c2: (abs(x[1] - x[2]) == 2) and ((x[1] != 2) or not (x[3] == 0));
s.t. c3: (x[3] == 0) or (abs(x[1] - x[2]) == 2);
s.t. c4: (x[1] != 2) ==> (x[1] != 2);
