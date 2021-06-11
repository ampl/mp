###  DATA  ###

set joints;
set bars within {i in joints, j in joints: i <> j};

				      # Definition of admissible structure:
				      # each bar connects two joints

param fixed symbolic in joints;       # Designated position of fixed support
param rolling symbolic in joints;     # Designated position of roller support

param density > 0;		      # Density of bar material
param yield_stress > 0;		      # Yield stress of bar material

param xpos {joints};		      # Horizontal positions of joints
param ypos {joints};		      # Vertical positions of joints

   check {(i,j) in bars}: xpos[i] <> xpos[j] or ypos[i] <> ypos[j];

param xload {joints};		      # Horizontal external loads on joints
param yload {joints};		      # Vertical external loads on joints

param length {(i,j) in bars} :=
		      sqrt ((xpos[j]-xpos[i])^2 + (ypos[j]-ypos[i])^2);

				      # Bar lengths calculated from positions

param xcos {(i,j) in bars} := (xpos[j]-xpos[i]) / length[i,j];
param ycos {(i,j) in bars} := (ypos[j]-ypos[i]) / length[i,j];

				      # Cosines of bar angles with
				      # horizontal and vertical axes

###  VARIABLES  ###

var Force {bars};		      # Forces on bars:
				      # positive in tension, negative in compression

###  OBJECTIVE  ###

minimize weight:  (density / yield_stress) *

   sum {(i,j) in bars} length[i,j] * <<0; -1,+1>> Force[i,j];

				      # Weight is proportional to length
				      # times absoluted value of force

###  CONSTRAINTS  ###

subject to xbal {k in joints: k <> fixed}:

     sum {(i,k) in bars} xcos[i,k] * Force[i,k]
   - sum {(k,j) in bars} xcos[k,j] * Force[k,j] = xload[k];

subject to ybal {k in joints: k <> fixed and k <> rolling}:

     sum {(i,k) in bars} ycos[i,k] * Force[i,k]
   - sum {(k,j) in bars} ycos[k,j] * Force[k,j] = yload[k];

				      # Net sum of forces must balance external
				      # load, horizontally and vertically
