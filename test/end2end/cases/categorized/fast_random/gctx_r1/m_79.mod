var x{1..5} integer >= -1 <= 1;

minimize o: -x[1] - x[2] - 4*x[3] + 3*x[4] - 2*x[5];

s.t. c1:
   (abs(x[5] - x[1]) <= 0)
   ==> ((abs(x[4] - x[5]) != 0) or (abs(x[5] - x[1]) <= 0));
s.t. c2: (if x[5] <= -1 then x[5] else x[3]) <= 1;
s.t. c3: (x[5] <= -1) or (x[5] <= -1);
