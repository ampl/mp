
Efficiency considerations
--------------------------

The goal of these extensions is to let you write models however you think about them, relying on the MP interface to convert them to the forms required by solvers. Nevertheless, there will be situations where one choice of formulation will lead to better solver performance than another. Here we collect some examples that are relevant to the current MP implementation. Future versions may automate some of these reformulations.


Simplification of logic
************************

.. code-block:: ampl

    var Flow {PRODUCTS,ARCS} >= 0;

    minimize TotalCost:
       sum {p in PRODUCTS, (i,j) in ARCS} var_cost[p,i,j] * Flow[p,i,j] +
       sum {(i,j) in ARCS}
          if exists {p in PRODUCTS} Flow[p,i,j] > 0 then fix_cost[i,j];

Each term ``Flow[p,i,j] > 0`` is converted separately, involving
a separate binary variable and implication constraint. But for a given
i and j, there exists a positive Flow[p,i,j] if and only if the sum of
all Flow[p,i,j] is positive. Thus an alternative formulation is given by:

.. code-block:: ampl

    minimize TotalCost:
       sum {p in PRODUCTS, (i,j) in ARCS} var_cost[p,i,j] * Flow[p,i,j] +
       sum {(i,j) in ARCS}
          if sum {p in PRODUCTS} Flow[p,i,j] > 0 then fix_cost[i,j];

By taking advantage of the solver's ability to work with linear expressions,
this form enables a substantially more concise reformulation.


Creation of common subexpressions
**********************************

.. code-block:: ampl

    subject to Shipment_Limits {(i,j) in ARCS}:
       sum {p in PRODUCTS} Flow[p,i,j] = 0 or
       min_ship <= sum {p in PRODUCTS} Flow[p,i,j] <= capacity[i,j];

    minimize TotalCost:
       sum {p in PRODUCTS, (i,j) in ARCS} var_cost[p,i,j] * Flow[p,i,j] +
       sum {(i,j) in ARCS}
          if sum {p in PRODUCTS} Flow[p,i,j] > 0 then fix_cost[i,j];

The constraint implies that if ``sum {p in PRODUCTS} Flow[p,i,j]`` is positive, then it must be at least equal to min_ship. Thus ``> 0`` can be replaced by ``>= min_ship`` in the objective expression:

.. code-block:: ampl

    minimize TotalCost:
       sum {p in PRODUCTS, (i,j) in ARCS} var_cost[p,i,j] * Flow[p,i,j] +
       sum {(i,j) in ARCS}
          if sum {p in PRODUCTS} Flow[p,i,j] >= min_ship then fix_cost[i,j];

As a result of this change, ``sum {p in PRODUCTS} Flow[p,i,j] >= min_ship`` is a subexpression in both the constraint and the objective, simplifying the converson of the model.


Bounds on subexpressions
*************************

.. code-block:: ampl

    var x {1..2} <= 2, >= -2;

    minimize Goldstein-Price:
       (1 + (x[1] + x[2] + 1)^2
          * (19 - 14*x[1] + 3*x[1]^2 - 14*x[2] + 6*x[1]*x[2] + 3*x[2]^2))
     * (30 + (2*x[1] - 3*x[2])^2
          * (18 - 32*x[1] + 12*x[1]^2 + 48*x[2] - 36*x[1]*x[2] + 27*x[2]^2));

Solver performance can often be improved by tightening bounds on the variables. In this example, the bounds on the variables also imply bounds on the subexpressions ``(x[1] + x[2] + 1)^2`` and ``(2*x[1] - 3*x[2])^2``. By defining auxiliary variables ``t1`` and ``t2`` equal to these subexpressions, their bounds can be communicated to the solver:

.. code-block:: ampl

    var t1 >= 0, <= 25;   subj to t1def: t1 = (x[1] + x[2] + 1)^2;
    var t2 >= 0, <= 100;  subj to t2def: t2 = (2*x[1] - 3*x[2])^2;

    minimize Goldstein-Price:
       (1 + t1
          * (19 - 14*x[1] + 3*x[1]^2 - 14*x[2] + 6*x[1]*x[2] + 3*x[2]^2))
     * (30 + t2
          * (18 - 32*x[1] + 12*x[1]^2 + 48*x[2] - 36*x[1]*x[2] + 27*x[2]^2));

These bounds are observed to substantially improve Gurobi's performance in this case.


