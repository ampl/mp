var x{1..4} integer >= -1 <= 2;

maximize o: -5*x[1] - 4*x[2] + 4*x[3] + 6*x[4];

s.t. c1:
   (x[4] + x[4] <= 0) and ((x[4] + x[4] <= 0) or not (max(x[1], x[3]) != 0));
s.t. c2: not ((x[4] + x[4] <= 0) and (max(x[1], x[3]) != 0));
s.t. c3:
   (max(x[1], x[3]) != 0) ==> ((max(x[1], x[3]) != 0) or (x[4] + x[4] <= 0));
s.t. c4:
   (max(x[1], x[3]) != 0) ==> ((x[4] + x[4] <= 0) or (max(x[1], x[3]) != 0));
