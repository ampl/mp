from amplpy import AMPL


def tryGetObjective(ampl):
    try: 
        return ampl.getCurrentObjective().value()
    except:
        try:
            # Get the first objective
            objs = ampl.getObjectives()
            obj = next(objs)[1]
            try:
                # Try as scalar objective
                return obj.value()
            except Exception as e: 
                pass
            for index, instance in obj:
                # Return the first instance if indexed
                return instance.value()
        except Exception as e:
            pass
    return None


ampl = AMPL()

ampl.reset()
ampl.eval("reset options;")


ampl.read("obj_suf_01__objpr_frac.mod")

ampl.solve(solver="gurobi", mp_options="multiobj=1")

print("Solve result =", ampl.getValue("solve_result"))
print("Solve time =", ampl.getValue("_solve_elapsed_time"))

print("Obj =", tryGetObjective(ampl))

# solver_time =ampl.get_value("Initial.time_solver")
# setup_time=ampl.get_value("Initial.time_setup")

print("Solve message:", ampl.getValue("solve_message"))



ampl.reset()
ampl.eval("reset options;")


ampl.read("obj_suf_01__objpr_frac.mod")

ampl.solve(solver="gurobi", mp_options="multiobj=1")

print("Solve message:", ampl.getValue("solve_message"))
print("Solve result =", ampl.solve_result)
