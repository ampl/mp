var x{1..6} integer >= 0 <= 1;

maximize o: 5*x[1] + 0*x[2] - 4*x[3] + x[4] - 3*x[5] - 3*x[6];

s.t. c1: (if x[2] != 0 then x[5] else x[5]) > 1;
s.t. c2: ((x[2] > 0) ==> (x[4] == 0)) and ((x[4] == 0) ==> (x[2] >= 0));
s.t. c3: max(x[6], x[6]) - min(x[6], x[3]) < 0;
s.t. c4: atleast 1 {i in 1..6} (x[i] >= 1);
