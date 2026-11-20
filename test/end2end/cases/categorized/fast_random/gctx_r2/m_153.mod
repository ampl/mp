var x{1..4} integer >= -1 <= 1;

minimize o: 5*x[1] - x[2] + 3*x[3] + 2*x[4];

s.t. c1: (x[2] > 0) ==> ((x[1] > 1) or (max(x[2], x[1]) < 0));
s.t. c2: (abs(x[2] - x[1]) < 1) <==> (abs(x[2] - x[1]) < 1);
