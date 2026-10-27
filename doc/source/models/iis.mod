var x >= 0;
var y >= 0;
var z >= 0;

subj to c1:
   x + y >= 3;

subj to c2:
   x + z <= 1;

subj to c3:
   y + z <= 1;

subj to c4:            # not needed for the infeasibility
   x + y + z <= 10;
