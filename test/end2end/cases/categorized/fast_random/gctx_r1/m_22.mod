var x{1..4} integer >= -1 <= 2;

minimize o:
   -2*x[1] + 6*x[2] - 2*x[3] - 2*x[4]
   + 0*(if abs(x[3] - x[1]) == 3 then 1 else 0);

s.t. c1: (x[1] == x[2]) and ((x[1] == x[2]) or not (x[1] == x[2]));
s.t. c2: (abs(x[3] - x[1]) == 3) and ((x[1] == x[2]) or not (x[1] == x[2]));
s.t. c3: not ((abs(x[3] - x[1]) == 3) and (x[1] == x[2]));
s.t. c4: not (abs(x[3] - x[1]) == 3);
