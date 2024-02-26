from scipy.sparse import csr_matrix
import numpy as np

from amplpy import modules
import nlwpy


class ModelBuilder:
    def __init__(self):
        self.prob_name_ = "nlwpy_prob"
        self.var_lb_ = [0, -3, 0, -1, -1, -2]
        self.var_ub_ = [0, 20, 1, 1e20, -1, 10]
        self.var_type_ = [0, 1, 1, 1, 0, 0]
        self.var_names_ = ["x1_4", "x2_6", "x3_5", "x4_3", "x5_1", "x6_2"]
        self.A_format_ = nlwpy.MatrixFormat.Rowwise
        self.A_ = np.array([[0, 1, 1, 1, 0, 1], [0, 1, -1, -1, 0, 1]])
        self.row_lb_ = [15, 10]
        self.row_ub_ = [15, np.inf]
        self.row_names_ = ["C1", "C2"]
        self.obj_sense_ = nlwpy.ObjSense.Minimize
        self.obj_c0_ = 3.24
        self.obj_c_ = [0, 1, 0, 0, 0, 0]
        self.Q_format_ = nlwpy.HessianFormat.Square
        self.Q_ = np.zeros([6, 6])
        self.Q_[3, 3] = 10
        self.Q_[3, 5] = 12
        self.Q_[4, 4] = 14
        self.obj_name_ = "obj[1]"

        ### Extra input
        self.ini_x_i_ = [0, 2]
        self.ini_x_v_ = [5, 4]
        self.ini_y_i_ = [0]
        self.ini_y_v_ = [-12]
        self.bas_x_ = [3, 4, 1, 4, 1, 3]  ### Basis statuses
        self.bas_y_ = [1, 1]

        ### Solution
        self.x_ref_ = [0, 5, 1, -1, -1, 10]
        self.obj_val_ref_ = -39.76

    def get_model(self):
        nlme = nlwpy.NLModel(self.prob_name_)

        nlme.SetCols(self.var_lb_, self.var_ub_, self.var_type_)
        nlme.SetColNames(self.var_names_)

        if self.A_ is not None:
            self.A_ = csr_matrix(self.A_)
            nlme.SetRows(
                self.row_lb_,
                self.row_ub_,
                self.A_format_,
                self.A_.indptr,
                self.A_.indices,
                self.A_.data,
            )
        nlme.SetRowNames(self.row_names_)

        nlme.SetLinearObjective(self.obj_sense_, self.obj_c0_, self.obj_c_)

        if self.Q_ is not None:
            self.Q_ = csr_matrix(self.Q_)
            nlme.SetHessian(
                self.Q_format_, self.Q_.indptr, self.Q_.indices, self.Q_.data
            )
        nlme.SetObjName(self.obj_name_)

        if len(self.ini_x_i_) > 0:
            nlme.SetWarmstart(self.ini_x_i_, self.ini_x_v_)
        if len(self.ini_y_i_) > 0:
            nlme.SetDualWarmstart(self.ini_y_i_, self.ini_y_v_)
        if len(self.bas_x_) > 0:
            suf = nlwpy.NLSuffix("status", 0, self.bas_x_)
            nlme.AddSuffix(suf)
        if len(self.bas_y_) > 0:
            suf = nlwpy.NLSuffix("status", 1, self.bas_y_)
            nlme.AddSuffix(suf)

        return nlme

    def check(self, sol):
        result = True
        if not self._approx_equal(sol.obj_val_, self.obj_val_ref_):
            print(
                "MIQP 1: wrong obj val ({:.17f} !~ {:.17f})".format(
                    sol.obj_val_, self.obj_val_ref_
                )
            )
            result = False

        for i in range(len(sol.x_)):
            if not self._approx_equal(self.x_ref_[i], sol.x_[i]):
                print(
                    "MIQP 1: wrong x[{}] ({:.17f} !~ {:.17f})".format(
                        i + 1, sol.x_[i], self.x_ref_[i]
                    )
                )
                result = False

        ### Printing suffixes.
        ### TODO replace by checking debug suffixes and ini guesses.
        suffixes = sol.suffixes_
        print("Number of suffixes returned:", len(suffixes))
        for suf in suffixes:
            print("    SUFFIX '{}' [{}]".format(suf.name_, suf.kind_))
            print("        Table:    ", suf.table_)
            print("        Values:   ", *suf.values_)
        sufMIPGapObj = suffixes.Find("relmipgap", 2)  ## should be 2+4
        if sufMIPGapObj is not None:
            print(
                "FOUND:    SUFFIX '{}' [{}].".format(
                    sufMIPGapObj.name_, sufMIPGapObj.kind_
                )
            )
        sufVarStatus = suffixes.Find("status", 0)
        if sufVarStatus is not None:
            print(
                "FOUND:    SUFFIX '{}' [{}].".format(
                    sufVarStatus.name_, sufVarStatus.kind_
                )
            )
        suffixes.clear()

        print(
            "MIQP 1: solution check {}, obj={:.17f}.".format(
                "ok" if result else "Failed", sol.obj_val_
            )
        )
        return result

    def _approx_equal(self, n, m):
        return abs(n - m) <= 1e-5 * min(1.0, abs(n) + abs(m))


def solve_and_check(solver, sopts, binary, stub):
    mb = ModelBuilder()
    nlme = mb.get_model()
    nlse = nlwpy.NLSolver()
    nlopts = nlwpy.MakeNLOptionsBasic_Default()
    assert 0 == nlopts.n_text_mode_
    assert 0 == nlopts.want_nl_comments_
    assert 1 == nlopts.flags_
    nlopts.n_text_mode_ = not binary
    nlopts.want_nl_comments_ = 1
    nlse.SetNLOptions(nlopts)
    nlse.SetFileStub(stub)
    sol = nlse.Solve(nlme, solver, sopts)
    if sol.solve_result_ > -2:  ## Some method for this?
        if not mb.check(sol):
            print("Solution check failed.")
            return False
    else:
        print(nlse.GetErrorMessage())
        return False

    return True


def find_solver():
    lst = modules.installed()
    assert "minos" in lst
    return modules.find("minos")


def test_binary_nl():
    solver = find_solver()
    sopts = ""
    binary = True
    stub = ""
    assert solve_and_check(solver, sopts, binary, stub)


def test_text_nl():
    solver = find_solver()
    sopts = ""
    binary = False
    stub = ""
    assert solve_and_check(solver, sopts, binary, stub)
