/**
 * impl6_if.mod.
 *
 * This example checks 2 things.
 * 1) that after appearing in CTX_POS,
 *    the context of x>=3.0 is updated to CTX_MIX
 *    (by the if-then in the objective).
 * 2) That the nonlinear discrete variable x' type
 *    is correctly recognized.
 */

var x integer, >=0.0, <=20.0;
var y >=0.0, <=25.0;
var b binary;

s.t. Con1: x+y <= 5.0;

## does not impose the case: s.t. Or1: x==0.0 or x>=3.0;

s.t. Impl1: b ==> x>=3.0;

minimize Semi: -x - b + if x>=3.0 then 10.33 else 0;

