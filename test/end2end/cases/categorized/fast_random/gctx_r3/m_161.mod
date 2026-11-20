var x{1..3} integer >= -1 <= 1;

minimize o:
   -3*x[1] - 3*x[2] - 2*x[3] + (if abs(x[3] - x[1]) != 0 then 1 else 0);

s.t. c1: not ((abs(x[3] - x[1]) != 0) and (x[1] == 0));
s.t. c2: (if abs(x[3] - x[1]) != 0 then x[2] else x[3]) < 0;
s.t. c3: (abs(x[3] - x[1]) != 0) <==> (x[1] == 0);
s.t. c4: (abs(x[3] - x[1]) != 0) or (x[1] == 0);
