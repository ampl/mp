var x{1..5} integer >= 0 <= 3;

maximize o: x[1] + 0*x[2] + x[3] + x[4] + 4*x[5];

s.t. c1: (x[1] > x[4]) ==> (x[2] < 3) else (x[1] > x[4]);
s.t. c2: (x[2] < 3) or (x[1] > x[4]);
s.t. c3: (x[2] < 3) <==> (x[2] < 3);
