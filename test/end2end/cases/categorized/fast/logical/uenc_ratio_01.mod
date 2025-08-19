var x in 1..9;
var y >= 1 <= 200;

minimize Obj: x+y;

s.t. Con: y>3 ==> (x==2 || x==6
       || x==5    ## With uenc:ratio>3, this should trigger uenc
);
