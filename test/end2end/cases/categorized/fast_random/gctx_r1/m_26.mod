var x{1..3} integer >= 0 <= 2;

maximize o: -2*x[1] + 4*x[2] + 3*x[3];

s.t. c1: (abs(x[3] - x[3]) != 2) <==> (x[1] == 2);
s.t. c2:
   (abs(x[3] - x[3]) != 2) and ((abs(x[3] - x[3]) != 2) or not (x[1] == 2));
s.t. c3: (if x[1] == 2 then x[2] else x[1]) < 1;
