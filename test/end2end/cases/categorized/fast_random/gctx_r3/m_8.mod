var x{1..5} integer >= -1 <= 2;

minimize o:
   3*x[1] + 5*x[2] - 4*x[3] + 0*x[4] + 3*x[5] + 3*(if x[5] >= 2 then 1 else 0);

s.t. c1: (x[5] >= 2) <==> (abs(x[5] - x[2]) > 2);
s.t. c2: (x[4] < x[1]) and ((x[4] < x[1]) or not (abs(x[5] - x[2]) > 2));
s.t. c3: (x[5] >= 2) and ((x[4] < x[1]) or not (abs(x[5] - x[2]) > 2));
s.t. c4: not ((abs(x[5] - x[2]) > 2) and (x[5] >= 2));
