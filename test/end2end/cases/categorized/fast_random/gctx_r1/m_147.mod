var x{1..3} integer >= -1 <= 2;

maximize o: -3*x[1] + 5*x[2] + 6*x[3];

s.t. c1: not (x[2] < 0);
s.t. c2: (x[2] + x[2] <= 2) and ((x[2] < 0) or not (x[1] == 2));
s.t. c3: (x[2] < 0) ==> (x[1] == 2);
s.t. c4: (x[2] < 0) <==> (x[1] == 2);
