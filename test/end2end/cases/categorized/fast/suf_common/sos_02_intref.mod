## Set integer values for .ref suffixes
## So they are passed as an integer suffix

var w integer >= 0 <= 3;
var x integer >= 0 <= 3;
var y integer >= 0 <= 3;
var z integer >= 0 <= 3;

maximize Obj: w + 2*x + 0.5*y + 3*z;
s.t. C1: w + x + y + z <= 5;

suffix sosno IN;
suffix ref IN;

let w.sosno := 1;
let x.sosno := 1;
let y.sosno := 1;
let z.sosno := 1;
let w.ref := 1;
let x.ref := 2;
let y.ref := 3;
let z.ref := 4;
