var x{1..6} integer >= -2 <= 0;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: not ((x[2] == -1) and (x[1] == 0));
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: (x[3] <= -1) <==> (x[2] + x[3] > -2);
s.t. c4: (if x[2] == -1 then x[4] else x[1]) > -2;
