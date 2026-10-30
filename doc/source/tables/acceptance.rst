.. This file is written by doc/scripts/acceptanceTable.py; do not edit.

.. dropdown:: Logical and conditional

   .. csv-table::
      :class: solver-support
      :header-rows: 1

      "Construct (option)","Gurobi 20261029","Xpress 20261029","CPLEX 20261029","COPT 20261029","MOSEK 20261029","HiGHS 20261029","SCIP 20261029","GCG 20261029","CBC 20250814","cuOpt 20261029","CP-SAT 20261029","BARONMP 20261029"
      "``and``, ``forall`` (``acc:and``)","C","C","–","–","–","–","C","–","–","–","C","–"
      "``or``, ``exists`` (``acc:or``)","C","C","–","–","–","–","C","–","–","–","C","–"
      "``not`` (``acc:not``)","–","–","–","–","–","–","–","–","–","–","–","–"
      "``==>`` ... ``else`` (``acc:impl``)","–","–","–","–","–","–","–","–","–","–","C","–"
      "``<==>`` (``acc:equiv``)","–","–","–","–","–","–","–","–","–","–","–","–"
      "``if`` ... ``then`` ... ``else`` (``acc:ifthen``)","–","–","–","–","–","–","–","–","–","–","C","–"
      "binary ``==>`` linear (``acc:indeq``)","C","C","C","C","C","–","C","–","–","–","C","–"
      "binary ``==>`` quadratic (``acc:indquadeq``)","–","–","–","–","–","–","–","–","–","–","–","–"

.. dropdown:: Counting and constraint programming

   .. csv-table::
      :class: solver-support
      :header-rows: 1

      "Construct (option)","Gurobi 20261029","Xpress 20261029","CPLEX 20261029","COPT 20261029","MOSEK 20261029","HiGHS 20261029","SCIP 20261029","GCG 20261029","CBC 20250814","cuOpt 20261029","CP-SAT 20261029","BARONMP 20261029"
      "``alldiff`` (``acc:alldiff``)","–","–","–","–","–","–","–","–","–","–","C","–"
      "``count``, ``atmost``, ``atleast``, ``exactly`` (``acc:count``)","–","–","–","–","–","–","–","–","–","–","C","–"
      "``numberof`` (``acc:numberofconst``)","–","–","–","–","–","–","–","–","–","–","C","–"

.. dropdown:: Piecewise-linear, min, max, abs, SOS

   .. csv-table::
      :class: solver-support
      :header-rows: 1

      "Construct (option)","Gurobi 20261029","Xpress 20261029","CPLEX 20261029","COPT 20261029","MOSEK 20261029","HiGHS 20261029","SCIP 20261029","GCG 20261029","CBC 20250814","cuOpt 20261029","CP-SAT 20261029","BARONMP 20261029"
      "``<<breakpoints; slopes>>`` (``acc:pl``)","C","C","C","–","–","–","–","–","–","–","–","–"
      "``max`` (``acc:max``)","C","CE","–","–","–","–","–","–","–","–","C","–"
      "``min`` (``acc:min``)","C","CE","–","–","–","–","–","–","–","–","C","–"
      "``abs`` (``acc:abs``)","C","CE","–","E","–","–","CE","–","–","–","C","–"
      "suffixes ``.sosno``, ``.ref`` (``acc:sos1``)","C","C","C","C","–","–","C","–","C","–","–","–"

.. dropdown:: Integer arithmetic and rounding

   .. csv-table::
      :class: solver-support
      :header-rows: 1

      "Construct (option)","Gurobi 20261029","Xpress 20261029","CPLEX 20261029","COPT 20261029","MOSEK 20261029","HiGHS 20261029","SCIP 20261029","GCG 20261029","CBC 20250814","cuOpt 20261029","CP-SAT 20261029","BARONMP 20261029"
      "``x*y`` (``acc:mul``)","–","–","–","–","–","–","–","–","","–","C","–"
      "``div`` (``acc:intdiv``)","–","–","–","–","–","–","–","–","","–","C","–"
      "``mod`` (``acc:mod``)","–","–","–","–","–","–","–","–","","–","C","–"
      "``floor`` (``acc:floor``)","–","–","–","–","–","–","–","–","","–","–","–"
      "``ceil`` (``acc:ceil``)","–","–","–","–","–","–","–","–","","–","–","–"
      "``round`` (``acc:round``)","–","–","–","–","–","–","–","–","","–","–","–"

.. dropdown:: Quadratic and conic

   .. csv-table::
      :class: solver-support
      :header-rows: 1

      "Construct (option)","Gurobi 20261029","Xpress 20261029","CPLEX 20261029","COPT 20261029","MOSEK 20261029","HiGHS 20261029","SCIP 20261029","GCG 20261029","CBC 20250814","cuOpt 20261029","CP-SAT 20261029","BARONMP 20261029"
      "quadratic ``<=``, ``=``, ``>=`` (``acc:quadle``)","C","C","C","C","C","–","C","–","–","–","–","C"
      "``lb <= quadratic <= ub`` (``acc:quadrange``)","–","–","–","–","C","–","C","–","–","–","–","C"
      "quadratic expression (``acc:quadfn``)","E","E","–","E","–","–","E","–","–","–","–","E"
      "second-order cone (``acc:quadcone``)","–","–","–","C","C","–","c","–","–","–","–","–"
      "rotated second-order cone (``acc:rotatedquadcone``)","–","–","–","C","C","–","–","–","–","–","–","–"
      "exponential cone (``acc:expcone``)","–","–","–","C","C","–","–","–","–","–","–","–"
      "power cone (``acc:powercone``)","–","–","–","–","–","–","–","–","–","–","–","–"
      "geometric cone (``acc:geomcone``)","–","–","–","–","–","–","–","–","–","–","–","–"

.. dropdown:: Nonlinear functions

   .. csv-table::
      :class: solver-support
      :header-rows: 1

      "Construct (option)","Gurobi 20261029","Xpress 20261029","CPLEX 20261029","COPT 20261029","MOSEK 20261029","HiGHS 20261029","SCIP 20261029","GCG 20261029","CBC 20250814","cuOpt 20261029","CP-SAT 20261029","BARONMP 20261029"
      "``exp`` (``acc:exp``)","CE","CE","–","E","–","–","CE","–","–","–","–","CE"
      "constant ``^`` expression (``acc:expa``)","CE","C","–","E","–","–","–","–","–","–","–","CE"
      "``log`` (``acc:log``)","CE","CE","–","E","–","–","CE","–","–","–","–","CE"
      "``log10`` (``acc:loga``)","CE","CE","–","E","–","–","CE","–","–","–","–","CE"
      "expression ``^`` constant (``acc:powc``)","CE","CE","–","E","–","–","CE","–","–","–","C","CE"
      "expression ``^`` expression (``acc:pow``)","–","–","–","–","–","–","–","–","–","–","–","–"
      "``/`` (``acc:div``)","E","CE","–","E","–","–","–","–","–","–","–","E"
      "``sin`` (``acc:sin``)","CE","CE","–","E","–","–","CE","–","–","–","–","–"
      "``cos`` (``acc:cos``)","CE","CE","–","E","–","–","CE","–","–","–","–","–"
      "``tan`` (``acc:tan``)","CE","CE","–","E","–","–","–","–","–","–","–","–"
      "``asin`` (``acc:asin``)","–","CE","–","E","–","–","–","–","–","–","–","–"
      "``acos`` (``acc:acos``)","–","CE","–","E","–","–","–","–","–","–","–","–"
      "``atan`` (``acc:atan``)","–","CE","–","E","–","–","–","–","–","–","–","–"
      "``sinh`` (``acc:sinh``)","–","c","–","E","–","–","–","–","–","–","–","–"
      "``cosh`` (``acc:cosh``)","–","c","–","E","–","–","–","–","–","–","–","–"
      "``tanh`` (``acc:tanh``)","E","c","–","E","–","–","–","–","–","–","–","–"
      "``asinh`` (``acc:asinh``)","–","c","–","E","–","–","–","–","–","–","–","–"
      "``acosh`` (``acc:acosh``)","–","c","–","E","–","–","–","–","–","–","–","–"
      "``atanh`` (``acc:atanh``)","–","c","–","E","–","–","–","–","–","–","–","–"
      "``1/(1+exp(-x))`` (``acc:logi``)","E","–","–","–","–","–","–","–","","–","–","–"
      "``sign(x)*abs(x)^a`` (``acc:signpowc``)","E","–","–","–","–","–","–","–","","–","–","–"
      "user functions (``acc:call``)","–","–","–","–","–","–","–","–","","–","–","–"

.. dropdown:: Complementarity

   .. csv-table::
      :class: solver-support
      :header-rows: 1

      "Construct (option)","Gurobi 20261029","Xpress 20261029","CPLEX 20261029","COPT 20261029","MOSEK 20261029","HiGHS 20261029","SCIP 20261029","GCG 20261029","CBC 20250814","cuOpt 20261029","CP-SAT 20261029","BARONMP 20261029"
      "``complements``, linear (``acc:compl``)","–","–","–","–","–","–","–","–","–","–","–","–"
      "``complements``, quadratic (``acc:complquad``)","–","–","–","–","–","–","–","–","–","–","–","–"
      "``complements``, nonlinear (``acc:nlcompl``)","–","–","–","–","–","–","–","–","–","–","–","–"

.. dropdown:: Linear and internal (technical)

   .. csv-table::
      :class: solver-support
      :header-rows: 1

      "Construct (option)","Gurobi 20261029","Xpress 20261029","CPLEX 20261029","COPT 20261029","MOSEK 20261029","HiGHS 20261029","SCIP 20261029","GCG 20261029","CBC 20250814","cuOpt 20261029","CP-SAT 20261029","BARONMP 20261029"
      "linear ``<=``, ``=``, ``>=`` (``acc:linle``)","C","C","C","C","C","C","C","C","C","C","C","C"
      "``lb <= linear <= ub`` (``acc:linrange``)","–","–","–","C","C","C","C","C","C","–","C","C"
      "linear expression (``acc:linfn``)","E","E","–","E","–","–","E","–","–","–","–","E"
      "conditional linear (``acc:condlineq``)","–","–","–","–","–","–","–","–","–","–","–","–"
      "conditional quadratic (``acc:condquadeq``)","–","–","–","–","–","–","–","–","–","–","–","–"
      "nonlinear constraint (``acc:nlcon``)","–","C","–","C","–","–","C","–","–","–","–","C"
      "nonlinear assignment (``acc:nlassigneq``)","C","C","–","C","–","–","C","–","–","–","–","C"
      "nonlinear logical (``acc:nllogcon``)","–","–","–","–","–","–","C","–","–","–","–","–"
      "nonlinear reified (``acc:nlreifequiv``)","–","–","–","–","–","–","C","–","–","–","–","–"
      "unary encoding (``acc:uenc``)","–","–","–","–","–","–","–","–","–","–","–","–"
      "SDP dot product (``acc:sdpdotprod``)","–","–","–","–","–","–","–","–","","–","–","–"
