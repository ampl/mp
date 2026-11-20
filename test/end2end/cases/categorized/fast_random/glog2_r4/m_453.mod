var x{1..4} integer >= 0 <= 3;

minimize o: -x[1] + 0*x[2] + 2*x[3] - x[4];

s.t. c1: ((x[1] != 3) ==> (x[1] == 0)) and ((x[1] == 0) ==> (x[1] >= 3));
s.t. c2: not ((x[1] >= 3) and (x[1] == 2));
s.t. c3: atleast 3 {i in 1..4} (x[i] >= 0);
s.t. c4: (x[4] != 3) or (x[2] > 2) or (x[1] == 3);
s.t. c5: x[3] == 0;
