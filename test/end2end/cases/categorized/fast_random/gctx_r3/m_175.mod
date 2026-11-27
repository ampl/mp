var x{1..5} integer >= 0 <= 3;

maximize o: -x[1] - 5*x[2] - 3*x[3] + x[4] + 0*x[5];

s.t. c1:
   (if max(x[1], x[3]) < 1 then x[3] else 0) + (if x[4] == 1 then x[5] else 0)
   <= 0;
s.t. c2: (max(x[1], x[3]) < 1) and ((x[4] == 1) or not (max(x[1], x[3]) < 1));
s.t. c3: (max(x[1], x[3]) < 1) <==> (max(x[1], x[3]) < 1);
s.t. c4: (max(x[1], x[3]) < 1) ==> ((max(x[1], x[3]) < 1) or (x[4] == 1));
