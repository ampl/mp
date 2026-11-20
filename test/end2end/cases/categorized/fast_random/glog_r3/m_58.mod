var x{1..3} integer >= -2 <= 2;

minimize o: -3*x[1] - 4*x[2] + 2*x[3];

s.t. c1: (if x[1] <= -2 then x[2] else x[1]) < -2;
s.t. c2: ((x[1] <= 2) ==> (x[2] > 0)) and ((x[2] > 0) ==> (x[1] >= 2));
s.t. c3: atleast 3 {i in 1..3} (x[i] >= 2);
