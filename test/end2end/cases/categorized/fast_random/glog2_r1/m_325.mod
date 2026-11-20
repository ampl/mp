var x{1..6} integer >= -2 <= 2;

maximize o: 5*x[1] + 4*x[2] + 4*x[3] + 2*x[4] - 4*x[5] + 2*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[5] != 0) ==> (x[6] <= -2) else (x[3] < -2);
s.t. c3: abs(x[2] - x[5]) != 0;
s.t. c4: abs(x[5] - x[5]) < 3;
s.t. c5: atleast 1 {i in 1..6} (x[i] >= -2);
s.t. c6: (x[2] != -2) <==> (x[2] + x[5] <= 1);
s.t. c7: (x[5] > 1) and (x[4] <= 1) ==> (x[5] >= 1);
