var x integer >= 5;
var y integer;

minimize TotalSum:
    x - 2*y;

subj to C1:
       -x + 21*y >=  2;

subj to C2:
     -3*x +  2*y <=  1;

subj to C3:
     20*x +    y <= 20;
