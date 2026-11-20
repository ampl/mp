var x{1..5} integer >= -1 <= 1;

maximize o: 0*x[1] - 5*x[2] + 2*x[3] - 2*x[4] + 6*x[5];

s.t. c1: not (x[4] + x[5] == 1);
s.t. c2:
   (abs(x[5] - x[2]) <= 1)
   ==> (abs(x[2] - x[5]) == 2)
   else (abs(x[2] - x[5]) == 2);
