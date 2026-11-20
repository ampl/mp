var x{1..5} integer >= 0 <= 3;

maximize o: 2*x[1] + 2*x[2] + 2*x[3] + 6*x[4] - 2*x[5];

s.t. c1: (abs(x[4] - x[5]) < 1) or (min(x[1], x[3]) > 2);
s.t. c2:
   (min(x[5], x[2]) == 0)
   and ((abs(x[4] - x[5]) < 1) or not (min(x[5], x[2]) == 0));
