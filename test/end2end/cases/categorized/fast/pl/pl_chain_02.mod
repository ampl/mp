# Chained PL terms with a start point, see pl_chain_01.mod.
# From this start point, the PL -> SOS2 -> MIP reformulation
# with end segments extended to +-1e6 gave wrong objectives.

var x1 := 0;   s.t. pin_x1: x1 = 0;
var x2 := -1;  s.t. pin_x2: x2 = -1;
var x3 := -1;  s.t. pin_x3: x3 = -1;

var y1 = <<-1,1; -2,0.5,0.5>> x1 - x2 - 3;             # -2
var y2 = <<-3,2; -2,0.5,-2>> y1 + x3;                  # -2

var t;
minimize o: t;                                         # 0
s.t. ct: t >= <<0,3; -1,-0.5,3>> y2 + y1;
