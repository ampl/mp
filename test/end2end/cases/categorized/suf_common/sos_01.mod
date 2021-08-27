
# -------------------------------------------------------------
# sos_01.mod
# -------------------------------------------------------------

param NVsos1 := 5;
param NVsos2 := 7;

var x {1..NVsos1} >= 0, <=3;
var y {1..NVsos2} >= 0, <=2.5;

maximize SomeElements:
    4*x[2] + x[5] + 4*y[3] + 4*y[2] + y[6];

subj to C1:
       sum {i in 1..NVsos1} x[i] + sum {i in 1..NVsos2} y[i] <= 30;
       
suffix sosno IN;
suffix ref IN;
       
let {i in 1..NVsos1} x[i].sosno := 1;   ## SOS1
let {i in 1..NVsos1} x[i].ref := (i+30)/2.0;

let {i in 1..NVsos2} y[i].sosno := -2;  ## SOS2
let {i in 1..NVsos2} y[i].ref := (i+30)/3.0;
