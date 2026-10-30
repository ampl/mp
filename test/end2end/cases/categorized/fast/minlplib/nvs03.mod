# MINLPLib instance nvs03, https://www.minlplib.org/nvs03.html
# (c) MINLPLib, licensed under CC BY 4.0,
# https://creativecommons.org/licenses/by/4.0/ ; only this header added.
#
#  MINLP written by GAMS Convert at 01/12/18 13:36:58
#  
#  Equation counts
#      Total        E        G        L        N        X        C        B
#          3        1        2        0        0        0        0        0
#  
#  Variable counts
#                   x        b        i      s1s      s2s       sc       si
#      Total     cont   binary  integer     sos1     sos2    scont     sint
#          3        1        0        2        0        0        0        0
#  FX      0        0        0        0        0        0        0        0
#  
#  Nonzero counts
#      Total    const       NL      DLL
#          7        4        3        0
# 
#  Reformulation has removed 1 variable and 1 equation


var i1 integer := 100, >= 0, <= 200;
var i2 integer := 100, >= 0, <= 200;

minimize obj: (-8 + i1)^2 + (-2 + i2)^2;

subject to

e1: -0.1*i1^2 + i2 >= 0;

e2:  - 0.333333333333333*i1 - i2 >= -4.5;
