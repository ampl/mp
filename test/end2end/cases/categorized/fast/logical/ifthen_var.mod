/**
  * Test expression map as well as if-then.
  * Requires cvt:cmp:eps > feastol
  * for AMPL to correctly compute obj value, see #102.
  */

var x >=-100, <= 200;
var y >=-300, <= 460;

var b: binary;

subj to ConImpl: b ==> 2*x + 3*y <= 5;

minimize TotalIf: if 2*x+3*y<=5 then b else -b-5;
