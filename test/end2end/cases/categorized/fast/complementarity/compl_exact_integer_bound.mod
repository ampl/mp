param B := 4503599627370497;
var x integer >= B-2, <= B+2;
var q >= 0;

minimize obj: x - B;

subject to KKT: x >= B complements q >= 0;
