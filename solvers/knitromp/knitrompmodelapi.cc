#include "knitrompmodelapi.h"


namespace mp {

void KnitrompModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {
  // TODO Allocate storage if needed:
  // auto n_linear_cons =
  //   flat_model_info->GetNumberOfConstraintsOfGroup(CG_LINEAR);
  // reallocate_linear_cons( n_linear_cons );
}

void KnitrompModelAPI::AddVariables(const VarArrayDef& v) {
    KN_add_vars(lp(), (int)v.size(), NULL);
    std::vector<int> xType(v.size());
    for (int i = 0; i < v.size(); ++i)
        xType[i] = v.ptype()[i] == var::INTEGER ? KN_VARTYPE_INTEGER : KN_VARTYPE_CONTINUOUS;
    KNITROMP_CCALL(KN_set_var_types_all(lp(), xType.data()));

	// Convert MP infinity to Knitro's representation of infinity
    auto convertBound = [](double bound) {
        if (bound == std::numeric_limits<double>::infinity())
            return KN_INFINITY;
        else if (bound == -std::numeric_limits<double>::infinity())
            return -KN_INFINITY;
        else
            return bound;
        };

    std::vector<double> bnds(v.size());
    std::transform(v.plb(), v.plb() + v.size(), bnds.begin(), convertBound);
    KNITROMP_CCALL(KN_set_var_lobnds_all(lp(), bnds.data()));

    std::transform(v.pub(), v.pub() + v.size(), bnds.begin(), convertBound);
    KNITROMP_CCALL(KN_set_var_upbnds_all(lp(), bnds.data()));

    std::vector<int> xLinear(v.size(), KN_VAR_LINEAR);
    KNITROMP_CCALL(KN_set_var_properties_all(lp(), xLinear.data()));
}


void KnitrompModelAPI::SetLinearObjective(int iobj, const LinearObjective& lo) {
    KNITROMP_CCALL(KN_add_obj_linear_struct(lp(), lo.num_terms(),
        lo.vars().data(), lo.coefs().data()));
	if (lo.obj_sense() == obj::MAX) {
        KNITROMP_CCALL(KN_set_obj_goal(lp(), KN_OBJGOAL_MAXIMIZE));
    }
}

void KnitrompModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
    KNITROMP_CCALL(KN_add_obj_linear_struct(lp(), qo.GetLinTerms().size(),
        qo.GetLinTerms().vars().data(), qo.GetLinTerms().coefs().data()));
    KNITROMP_CCALL(KN_add_obj_quadratic_struct(lp(), qo.GetQPTerms().size(),
        qo.GetQPTerms().vars1().data(), qo.GetQPTerms().vars2().data(),
		qo.GetQPTerms().coefs().data()));

    if (qo.obj_sense() == obj::MAX) {
        KNITROMP_CCALL(KN_set_obj_goal(lp(), KN_OBJGOAL_MAXIMIZE));
    }
}



void KnitrompModelAPI::AddConstraint(const LinConLE& lc) {
    int index;
    KN_add_con(lp(), &index);
    KN_add_con_linear_struct_one(lp(), lc.size(), index, lc.pvars(), lc.pcoefs());
    KN_set_con_upbnd(lp(), index, lc.rhs());
}
void KnitrompModelAPI::AddConstraint(const LinConEQ& lc) {
    int index;
    KN_add_con(lp(), &index);
    KN_add_con_linear_struct_one(lp(), lc.size(), index, lc.pvars(), lc.pcoefs());
    KN_set_con_eqbnd(lp(), index, lc.rhs());
}
void KnitrompModelAPI::AddConstraint(const LinConGE& lc) {
    int index;
    KN_add_con(lp(), &index);
    KN_add_con_linear_struct_one(lp(), lc.size(), index, lc.pvars(), lc.pcoefs());
    KN_set_con_lobnd(lp(), index, lc.rhs());
}

void KnitrompModelAPI::AddConstraint( const QuadConLE& qc ) {
    int index;
    KN_add_con(lp(), &index);
    auto lc = qc.GetLinTerms();
    KN_add_con_linear_struct_one(lp(), lc.size(), index, lc.pvars(), lc.pcoefs());
    KN_add_con_quadratic_struct_one(lp(), qc.GetQPTerms().size(), index, qc.GetQPTerms().vars1().data(),
        qc.GetQPTerms().vars2().data(), qc.GetQPTerms().coefs().data());
    KN_set_con_upbnd(lp(), index, qc.rhs());
}

void KnitrompModelAPI::AddConstraint( const QuadConEQ& qc ) {
    int index;
    KN_add_con(lp(), &index);
    auto lc = qc.GetLinTerms();
    KN_add_con_linear_struct_one(lp(), lc.size(), index, lc.pvars(), lc.pcoefs());
    KN_add_con_quadratic_struct_one(lp(), qc.GetQPTerms().size(), index, qc.GetQPTerms().vars1().data(),
		qc.GetQPTerms().vars2().data(), qc.GetQPTerms().coefs().data());
    KN_set_con_eqbnd(lp(), index, qc.rhs());
}

void KnitrompModelAPI::AddConstraint( const QuadConGE& qc ) {
    int index;
    KN_add_con(lp(), &index);
    auto lc = qc.GetLinTerms();
    KN_add_con_linear_struct_one(lp(), lc.size(), index, lc.pvars(), lc.pcoefs());
    KN_add_con_quadratic_struct_one(lp(), qc.GetQPTerms().size(), index, qc.GetQPTerms().vars1().data(),
        qc.GetQPTerms().vars2().data(), qc.GetQPTerms().coefs().data());
    KN_set_con_lobnd(lp(), index, qc.rhs());

}

template<typename UnaryOp>
void KnitrompModelAPI::AddUnaryNonLinearConstraint(int resultVar, int argVar, UnaryOp op) {
    std::vector<CppAD::AD<double>> X(1);
    CppAD::Independent(X);
    std::vector<CppAD::AD<double>> Y(1);
    Y[0] = op(X[0]);  
    CppAD::ADFun<double> tape(X, Y);
    AddNonlinearConstraint(resultVar, { argVar }, std::move(tape));
}

// Trigonometric and hyperbolic functions
void KnitrompModelAPI::AddConstraint(const SinConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::sin(x); });
}
void KnitrompModelAPI::AddConstraint(const CosConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::cos(x); });
}
void KnitrompModelAPI::AddConstraint(const TanConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::tan(x); });
}
void KnitrompModelAPI::AddConstraint(const AsinConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::asin(x); });
}
void KnitrompModelAPI::AddConstraint(const AcosConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::acos(x); });
}
void KnitrompModelAPI::AddConstraint(const AtanConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::atan(x); });
}
void KnitrompModelAPI::AddConstraint(const SinhConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::sinh(x); });
}
void KnitrompModelAPI::AddConstraint(const CoshConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::cosh(x); });
}
void KnitrompModelAPI::AddConstraint(const TanhConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::tanh(x); });
}
void KnitrompModelAPI::AddConstraint(const AsinhConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::asinh(x); });
}
void KnitrompModelAPI::AddConstraint(const AcoshConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::acosh(x); });
}
void KnitrompModelAPI::AddConstraint(const AtanhConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::atanh(x); });
}

void KnitrompModelAPI::AddConstraint(const LogConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::log(x); });
}
void KnitrompModelAPI::AddConstraint(const ExpConstraint& c) {
    AddUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0],
        [](const CppAD::AD<double>& x) { return CppAD::exp(x); });
}

void KnitrompModelAPI::AddConstraint(const LogAConstraint& c) {
    std::vector<CppAD::AD<double>> X(1);
    X[0] = 1.0;
    CppAD::Independent(X);
    std::vector<CppAD::AD<double>> Y(1);
    Y[0] = CppAD::log(X[0]) / log(c.GetParameters()[0]);
    CppAD::ADFun<double> tape(X, Y);
    std::vector<int> args = { c.GetArguments()[0] };
    AddNonlinearConstraint(c.GetResultVar(), args, std::move(tape));
}
void KnitrompModelAPI::AddConstraint(const PowConstraint& c) {
    std::vector<CppAD::AD<double>> X(2);
    X[0] = 1.0;
    X[1] = 1.0;
    CppAD::Independent(X);
    std::vector<CppAD::AD<double>> Y(1);
    Y[0] = CppAD::pow(X[0], X[1]);
    CppAD::ADFun<double> tape(X, Y);
    std::vector<int> args = { c.GetArguments()[0], c.GetArguments()[1] };
    AddNonlinearConstraint(c.GetResultVar(), args, std::move(tape));
}
void KnitrompModelAPI::AddConstraint(const DivConstraint& c) {
	std::vector<CppAD::AD<double>> X(2);
    // Initialize with dummy values for tape creation, not sure why
    // it matters here and not in unary constraints
	X[0] = 1.0;  
	X[1] = 1.0;
	CppAD::Independent(X);
	std::vector<CppAD::AD<double>> Y(1);
	Y[0] = X[0] / X[1];  
	CppAD::ADFun<double> tape(X, Y);
	std::vector<int> args = { c.GetArguments()[0], c.GetArguments()[1] };
	AddNonlinearConstraint(c.GetResultVar(), args, std::move(tape));
}

/// Helper to add an explicited linear constraint  
/// Stores the constraint data in nlConstraints_ for later evaluation in the callback
void KnitrompModelAPI::AddNonlinearConstraint(
    int resultVar,
    const std::vector<int>& argVars,
    CppAD::ADFun<double>&& tape) {

    // Add a nonlinear constraint to Knitro: resultVar - f(argVars) = 0
    int conIndex;
    KNITROMP_CCALL(KN_add_con(lp(), &conIndex));

    // Set as equality constraint with bound 0
    KNITROMP_CCALL(KN_set_con_eqbnd(lp(), conIndex, 0.0));

    // nonlinear - resultvar = 0
    int indexVar = resultVar;
    double coef = -1.0;
    KNITROMP_CCALL(KN_add_con_linear_struct_one(
        lp(), 1, conIndex, &indexVar, &coef));

    // Store the nonlinear constraint data
    NonlinearConstraintData data;
    data.tape = std::move(tape);
    data.tape.check_for_nan(false);
    data.resultVar = resultVar;
    data.argVars = argVars;
    data.knitroConIndex = conIndex;
    int nlConIndex = static_cast<int>(nlConstraints_.size());
    nlConstraints_.push_back(std::move(data));

    // Map Knitro constraint index to our internal index
    knitroConIndexToNLConIndex_[conIndex] = nlConIndex;
}

int KnitrompModelAPI::evalNonlinearConstraint(
	KN_context_ptr kc,
	CB_context_ptr cb,
	KN_eval_request_ptr const evalRequest,
	KN_eval_result_ptr const evalResult,
	void* const userParams) {

	auto* self = static_cast<KnitrompModelAPI*>(userParams);
	const double* x = evalRequest->x;

	if (evalRequest->type == KN_RC_EVALFC) {
		// Evaluate all nonlinear constraints
		for (size_t i = 0; i < self->nlConstraints_.size(); ++i) {
			auto& nlCon = self->nlConstraints_[i];

			// Get argument values
			std::vector<double> args(nlCon.argVars.size());
			for (size_t j = 0; j < nlCon.argVars.size(); ++j) {
				args[j] = x[nlCon.argVars[j]];
			}

            double banana;
            if (args[0] != 0) {
                banana = log(args[0]);
                printf("Point = %f, value = %f ", args[0], banana);
            }

			std::vector<double> y = nlCon.tape.Forward(0, args);
			printf("Value = %f\n", y[0]);
            if (std::isnan(y[0]) || std::isinf(y[0])) {
                evalResult->c[i] = KN_RC_EVAL_ERR;
            }
            else {
                evalResult->c[i] = y[0];
            }
		}
	}
	else if (evalRequest->type == KN_RC_EVALGA) {
		// Evaluate Jacobian for all nonlinear constraints
		for (size_t i = 0; i < self->nlConstraints_.size(); ++i) {
			auto& nlCon = self->nlConstraints_[i];

			// Get argument values
			std::vector<double> args(nlCon.argVars.size());
			for (size_t j = 0; j < nlCon.argVars.size(); ++j) {
				args[j] = x[nlCon.argVars[j]];
			}

			// Compute Jacobian using CppAD
          

			std::vector<double> jac = nlCon.tape.Jacobian(args);
	
			// Store Jacobian elements
            for (size_t j = 0; j < jac.size(); ++j) {
                if (std::isnan(jac[j]) || std::isinf(jac[j])) {
                    evalResult->jac[j] = KN_RC_EVAL_ERR;
                }
                else {
                    evalResult->jac[j] = jac[j];
                }
            }

		}
	}

	return KN_RC_EVALFC;
}

void KnitrompModelAPI::FinishProblemModificationPhase() {

    // Register callback for nonlinear constraints
    if (!nlConstraints_.empty()) {
        // Collect all constraint indices
        std::vector<int> conIndices;
        std::vector<int> jacIndexCons;
        std::vector<int> jacIndexVars;

        for (const auto& nlCon : nlConstraints_) {
            conIndices.push_back(nlCon.knitroConIndex);
            // Add Jacobian structure for this constraint
            for (int argVar : nlCon.argVars) {
                jacIndexCons.push_back(nlCon.knitroConIndex);
                jacIndexVars.push_back(argVar);
            }
        }

        CB_context_ptr cbContext;
        KNITROMP_CCALL(KN_add_eval_callback(
            lp(), 
            KNFALSE,  // evalObj - not evaluating objective
            static_cast<int>(conIndices.size()),
            conIndices.data(),
            evalNonlinearConstraint,
            &cbContext));

        // Store pointer to this object as user data
        KNITROMP_CCALL(KN_set_cb_user_params(
            lp(), cbContext, static_cast<void*>(this)));

        // Set up gradient structure for all constraints
        KNITROMP_CCALL(KN_set_cb_grad(
            lp(), cbContext,
            0, nullptr,  // No objective gradient
            static_cast<int>(jacIndexVars.size()),
            jacIndexCons.data(),
            jacIndexVars.data(),
            evalNonlinearConstraint));
    }
}


} // namespace mp
