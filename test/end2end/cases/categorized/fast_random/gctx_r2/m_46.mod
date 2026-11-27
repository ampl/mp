var x{1..5} integer >= 0 <= 3;

maximize o: -5*x[1] + 6*x[2] + 0*x[3] + x[4] + 4*x[5];

s.t. c1: (x[2] + x[3] < 2) ==> (x[4] < 1) else (x[4] < 1);
s.t. c2:
   (abs(x[3] - x[1]) < 1) and ((x[2] + x[3] < 2) or not (abs(x[3] - x[1]) < 1));
s.t. c3: not ((abs(x[3] - x[1]) < 1) and (x[4] < 1));
s.t. c4: (x[2] + x[3] < 2) <==> (x[4] < 1);
