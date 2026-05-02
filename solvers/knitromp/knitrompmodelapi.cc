#include "knitrompmodelapi.h"


namespace mp {

    void KnitrompModelAPI::InitProblemModificationPhase(
        const FlatModelInfo*) {
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

        if (v.pnames())
            KNITROMP_CCALL(KN_set_var_names_all(lp(), v.pnames()));
    }


    void KnitrompModelAPI::SetLinearObjective(int iobj, const LinearObjective& lo) {
        KNITROMP_CCALL(KN_add_obj_linear_struct(lp(), lo.num_terms(),
            lo.vars().data(), lo.coefs().data()));
        if (lo.obj_sense() == obj::MAX) {
            KNITROMP_CCALL(KN_set_obj_goal(lp(), KN_OBJGOAL_MAXIMIZE));
        }
        if (lo.name())
            KNITROMP_CCALL(KN_set_obj_name(lp(), lo.name()));

        if (printProblem)
        {
            PrintConstraint p(0, 0, lo.vars(), lo.coefs());
            printConstraints_[-1] = p;
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

        nonlinearVars_.insert(qo.GetQPTerms().vars1().begin(), qo.GetQPTerms().vars1().end());
        nonlinearVars_.insert(qo.GetQPTerms().vars2().begin(), qo.GetQPTerms().vars2().end());

        if (qo.name())
            KNITROMP_CCALL(KN_set_obj_name(lp(), qo.name()));

        if (printProblem)
        {
            auto lo = qo.GetLinTerms();
            auto qt = qo.GetQPTerms();
            PrintConstraint p(0, 0, lo.vars(), lo.coefs(), qt.vars1(), qt.vars2(), qt.coefs());
            printConstraints_[-1] = p;
        }
    }


    // Helper function for adding linear constraints
    template<typename LinConType>
    void KnitrompModelAPI::AddLinearConstraintHelper(
        const LinConType& lc,
        double lb,
        double ub) {

        int index;
        KN_add_con(lp(), &index);
        KN_add_con_linear_struct_one(lp(), lc.size(), index, lc.pvars(), lc.pcoefs());

        // Set bounds based on constraint type
        if (lb == ub) {
            // Equality constraint
            KN_set_con_eqbnd(lp(), index, lb);
        }
        else {
            // Inequality constraint(s)
            if (lb > -KN_INFINITY)
                KN_set_con_lobnd(lp(), index, lb);
            if (ub < KN_INFINITY)
                KN_set_con_upbnd(lp(), index, ub);
        }

        // Set constraint name if provided
        if (lc.name())
            KNITROMP_CCALL(KN_set_con_name(lp(), index, lc.name()));

        // Store for printing if needed
        if (printProblem) {
            PrintConstraint p(lb, ub, lc.vars(), lc.coefs());
            printConstraints_[NumCons() - 1] = p;
        }
    }



    void KnitrompModelAPI::AddConstraint(const LinConLE& lc) {
        AddLinearConstraintHelper(lc, -KN_INFINITY, lc.rhs());
    }

    void KnitrompModelAPI::AddConstraint(const LinConEQ& lc) {
        AddLinearConstraintHelper(lc, lc.rhs(), lc.rhs());
    }

    void KnitrompModelAPI::AddConstraint(const LinConGE& lc) {
        AddLinearConstraintHelper(lc, lc.rhs(), KN_INFINITY);
    }

    // Helper function for adding quadratic constraints
    template<typename QuadConType>
    void KnitrompModelAPI::AddQuadraticConstraintHelper(
        const QuadConType& qc,
        double lb,
        double ub) {

        int index;
        KN_add_con(lp(), &index);

        // Add linear terms
        auto lc = qc.GetLinTerms();
        KN_add_con_linear_struct_one(lp(), lc.size(), index, lc.pvars(), lc.pcoefs());

        // Add quadratic terms
        const auto& qt = qc.GetQPTerms();
        KN_add_con_quadratic_struct_one(lp(), qt.size(), index,
            qt.vars1().data(), qt.vars2().data(), qt.coefs().data());

        // Set bounds based on constraint type
        if (lb == ub) {
            KN_set_con_eqbnd(lp(), index, lb);
        }
        else {
            if (lb > -KN_INFINITY)
                KN_set_con_lobnd(lp(), index, lb);
            if (ub < KN_INFINITY)
                KN_set_con_upbnd(lp(), index, ub);
        }

        // Mark variables appearing in quadratic terms as nonlinear
        // Use index-based insertion to avoid iterator issues with ArrayRef
        for (size_t i = 0; i < qt.size(); ++i) {
            nonlinearVars_.insert(qt.vars1()[i]);
            nonlinearVars_.insert(qt.vars2()[i]);
        }

        // Set constraint name if provided
        if (qc.name())
            KNITROMP_CCALL(KN_set_con_name(lp(), index, qc.name()));

        // Store for printing if needed
        if (printProblem) {
            PrintConstraint p(lb, ub, lc.vars(), lc.coefs(),
                qt.vars1(), qt.vars2(), qt.coefs());
            printConstraints_[NumCons() - 1] = p;
        }
    }

    void KnitrompModelAPI::AddConstraint(const QuadConLE& qc) {
        AddQuadraticConstraintHelper(qc, -KN_INFINITY, qc.rhs());
    }

    void KnitrompModelAPI::AddConstraint(const QuadConEQ& qc) {
        AddQuadraticConstraintHelper(qc, qc.rhs(), qc.rhs());
    }

    void KnitrompModelAPI::AddConstraint(const QuadConGE& qc) {
        AddQuadraticConstraintHelper(qc, qc.rhs(), KN_INFINITY);
    }

    
    void KnitrompModelAPI::addUnaryNonLinearConstraint(int resultVar, int argVar, ExpressionData::OpType op, const char* name) {
		auto exp = ExpressionData::MakeUnaryExpr(op, ExpressionData::MakeVarExpr(argVar));
        addAssignNonlinearConstraint(EQ, 0.0, resultVar, exp, name);
    }

    // Trigonometric and hyperbolic functions       
    void KnitrompModelAPI::AddConstraint(const SinConstraint& c) {
		addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::SIN, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const CosConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::COS, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const TanConstraint& c) {
		addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::TAN, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const AsinConstraint& c) {
		addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::ASIN, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const AcosConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::ACOS, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const AtanConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::ATAN, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const SinhConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::SINH, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const CoshConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::COSH, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const TanhConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::TANH, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const AsinhConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::ASINH, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const AcoshConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::ACOSH, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const AtanhConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::ATANH, c.name());
    }

    void KnitrompModelAPI::AddConstraint(const LogConstraint& c) {
        addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::LOG, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const ExpConstraint& c) {
		addUnaryNonLinearConstraint(c.GetResultVar(), c.GetArguments()[0], ExpressionData::EXP, c.name());
    }

    void KnitrompModelAPI::AddConstraint(const LogAConstraint& c) {
        auto log_x = ExpressionData::MakeUnaryExpr(ExpressionData::LOG, ExpressionData::MakeVarExpr(c.GetArguments()[0]));
        double base = c.GetParameters()[0];
        double log_base = std::log(base);
        auto log_base_expr = ExpressionData::MakeConstantExpr(log_base);
        auto exp = ExpressionData::BinaryOp(ExpressionData::DIV, log_x, log_base_expr);
        addAssignNonlinearConstraint(EQ, 0.0, c.GetResultVar(), exp, c.name());
    }
    void KnitrompModelAPI::AddConstraint(const PowConstraint& c) {
        auto exp = ExpressionData::BinaryOp(ExpressionData::POW,
            ExpressionData::MakeVarExpr(c.GetArguments()[0]),
            ExpressionData::MakeVarExpr(c.GetArguments()[1]));
        addAssignNonlinearConstraint(EQ, 0.0, c.GetResultVar(), exp, c.name());


    }
    void KnitrompModelAPI::AddConstraint(const DivConstraint& c) {
        auto exp =ExpressionData::BinaryOp(ExpressionData::DIV, 
            ExpressionData::MakeVarExpr(c.GetArguments()[0]),
            ExpressionData::MakeVarExpr(c.GetArguments()[1]));
		addAssignNonlinearConstraint(EQ, 0.0, c.GetResultVar(), exp, c.name());
    }


    // Store the nonlinear constraint data for evaluation.
    // To be called after adding the constraint to Knitro.
    void  KnitrompModelAPI::storeNonLinearData(int conIndex, const ExpressionData& exp,
        int resultVar, const char* name)
    {

        NonlinearConstraintData data;
        auto [tape, argVars] = createTape(exp);
        data.tape = std::make_shared<CppAD::ADFun<double>>(std::move(tape));
        data.tape->check_for_nan(false);
        data.resultVar = resultVar;
        data.argVars = argVars;
        data.knitroConIndex = conIndex;
        if (printProblem)
            data.originalExpr = std::make_shared<ExpressionData>(exp);
        
        if (useHessian) {
            // Compute hessian sparsity manually
			// TODO use CPPAD's ForSparseHes, although in my first tests it seems to be 
            // giving incorrect results (missing some entries in the sparsity pattern)
            auto hs = ExpressionData::computeHessianSparsityInfo(exp);

            data.hessianSparsity.clear();

            for (const auto& [globalRow, globalCol] : hs.hessPairs) {
                auto rowIt = std::find(argVars.begin(), argVars.end(), globalRow);
                auto colIt = std::find(argVars.begin(), argVars.end(), globalCol);

                if (rowIt == argVars.end() || colIt == argVars.end())
                    continue;

                int localRow = static_cast<int>(std::distance(argVars.begin(), rowIt));
                int localCol = static_cast<int>(std::distance(argVars.begin(), colIt));

                if (localRow < localCol)
                    std::swap(localRow, localCol);

                data.hessianSparsity.push_back({ localRow, localCol });
            }
        }
        if (conIndex == -1) { // objective 
            nlObjective_ = data;
        }
        else {
            nlConstraints_.push_back(std::move(data));
            nonLinearConstraintsIndex_.push_back(conIndex);
            knitroToNLConstraintIndex_[conIndex] = nlConstraints_.size() - 1;
        }
        if (name != nullptr) {
            if (conIndex == -1)
                KNITROMP_CCALL(KN_set_obj_name(lp(), name));
            else
                KNITROMP_CCALL(KN_set_con_name(lp(), conIndex, name));

        }
    }



	/// Overload from expression class, builds the tape from the expression
    // and uses AddAssignNonLinearConstraint
    void KnitrompModelAPI::addAssignNonlinearConstraint(
        ConstraintType type, double rhs,
        int resultVar, const ExpressionData &exp, const char* name = nullptr) {
        
        // Add a nonlinear constraint to Knitro: f(argVars) - resultVar <type> rhs
        int conIndex;
        KNITROMP_CCALL(KN_add_con(lp(), &conIndex));
        if (type == ConstraintType::LE) {
            KNITROMP_CCALL(KN_set_con_upbnd(lp(), conIndex, rhs));
        }
        else if (type == ConstraintType::GE) {
            KNITROMP_CCALL(KN_set_con_lobnd(lp(), conIndex, rhs));
        }
        else
        {
            // Set as equality constraint with bound 0
            KNITROMP_CCALL(KN_set_con_eqbnd(lp(), conIndex, rhs));
        }

        int indexVar = resultVar;
        double coef = -1.0;
        KNITROMP_CCALL(KN_add_con_linear_struct_one(
            lp(), 1, conIndex, &indexVar, &coef));
        storeNonLinearData(conIndex, exp, resultVar, name);
    }

    int KnitrompModelAPI::evalNonlinearConstraint(
        KN_context_ptr kc,
        CB_context_ptr cb,
        KN_eval_request_ptr const evalRequest,
        KN_eval_result_ptr const evalResult,
        void* const userParams) {

        auto* self = static_cast<KnitrompModelAPI*>(userParams);
        const double* x = evalRequest->x;

        if ((evalRequest->type == KN_RC_EVALFC)  || (evalRequest->type == KN_RC_EVALGA) || (evalRequest->type == KN_RC_EVALFCGA)) {

            // Evaluate all nonlinear constraints
            for (size_t i = 0; i < self->nlConstraints_.size(); ++i) {
                auto& nlCon = self->nlConstraints_[i];

                // Get argument values
                std::vector<double> args(nlCon.argVars.size());
                for (size_t j = 0; j < nlCon.argVars.size(); ++j) {
                    args[j] = x[nlCon.argVars[j]];
                }
                // Want function value?
                if ((evalRequest->type == KN_RC_EVALFC) || (evalRequest->type == KN_RC_EVALFCGA))
                {
                    std::vector<double> y = nlCon.tape->Forward(0, args);
                    if (std::isnan(y[0]) || std::isinf(y[0])) {
                        evalResult->c[i] = KN_RC_EVAL_ERR;
                    }
                    else {
                        evalResult->c[i] = y[0];
                    }
                }

                // Want gradient/Jacobian?
                if ((evalRequest->type == KN_RC_EVALGA) || (evalRequest->type == KN_RC_EVALFCGA)) {
                    std::vector<double> jac = nlCon.tape->Jacobian(args);
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

			// Evaluate objective if needed
            if (self->nlObjective_.has_value()) {
                // Get argument values
                std::vector<double> args(self->nlObjective_->argVars.size());
                for (size_t j = 0; j < self->nlObjective_->argVars.size(); ++j) {
                    args[j] = x[self->nlObjective_->argVars[j]];
                }
                // Want function value?
                if ((evalRequest->type == KN_RC_EVALFC) || (evalRequest->type == KN_RC_EVALFCGA)) {
                    std::vector<double> y = self->nlObjective_->tape->Forward(0, args);
                    if (std::isnan(y[0]) || std::isinf(y[0])) {
                        *evalResult->obj = KN_RC_EVAL_ERR;
                    }
                    else {
                        *evalResult->obj = y[0];
                    }
                }

                // Want gradient?
                if ((evalRequest->type == KN_RC_EVALGA) || (evalRequest->type == KN_RC_EVALFCGA)) {
                    // Compute Jacobian using CppAD
                    std::vector<double> grad = self->nlObjective_->tape->Jacobian(args);
                    // Store Jacobian elements
                    for (size_t j = 0; j < grad.size(); ++j) {
                        if (std::isnan(grad[j]) || std::isinf(grad[j]))
                            evalResult->objGrad[j] = KN_RC_EVAL_ERR;
                        else
                            evalResult->objGrad[j] = grad[j];
                    }

                }
            }
        }
        return 0;
    }

    int KnitrompModelAPI::evalHessian(
        KN_context_ptr kc,
        CB_context_ptr cb,
        KN_eval_request_ptr const evalRequest,
        KN_eval_result_ptr const evalResult,
        void* const userParams) {

        auto* self = static_cast<KnitrompModelAPI*>(userParams);
        const double* x = evalRequest->x;
        const double* lambda = evalRequest->lambda;
        const double sigma = *(evalRequest->sigma);
        double* hess = evalResult->hess;

        if (evalRequest->type == KN_RC_EVALH) {
            // Initialize to zero - we'll accumulate contributions
            std::fill(hess, hess + self->hessSparsityMap_.size(), 0.0);

            // Helper to accumulate Hessian contributions
            auto accumulateHessian = [&hess](
                const NonlinearConstraintData& nlData,
                const double* x,
                double weight) {

                    // Get argument values for this tape
                    std::vector<double> args(nlData.argVars.size());
                    for (size_t j = 0; j < nlData.argVars.size(); ++j) {
                        args[j] = x[nlData.argVars[j]];
                    }

                    // Compute Hessian using CppAD, weighted
                    std::vector<double> w(1);
                    w[0] = weight;
                    std::vector<double> hesLocal = nlData.tape->Hessian(args, w);

                    // Add local Hessian elements to global array
                    for (size_t k = 0; k < nlData.hessianSparsity.size(); ++k) {
                        auto [row, col] = nlData.hessianSparsity[k];

                        // Convert local sparsity coordinates to local Hessian array index
                        // CppAD returns lower triangle in column-major (packed) format
                        size_t localIdx = row * (row + 1) / 2 + col;

                        // Get global Hessian array index
                        size_t globalIdx = nlData.hessianGlobalIndices[k];

                        // Accumulate (handles duplicates from different constraints/obj)
                        if (!std::isnan(hesLocal[localIdx]) && !std::isinf(hesLocal[localIdx])) {
                            hess[globalIdx] += hesLocal[localIdx];
                        }
                    }
                };

            // Add objective Hessian contribution
            if (self->nlObjective_.has_value()) {
                accumulateHessian(*self->nlObjective_, x, sigma);
            }

            // Add constraint Hessian contributions
            for (const auto& nlCon : self->nlConstraints_) {
                accumulateHessian(nlCon, x, lambda[nlCon.knitroConIndex]);
            }

            return 0;
        }
        return KN_RC_CALLBACK_ERR;
    }

    void KnitrompModelAPI::markLinearVariables() {
		int numVars = NumVars();

        for (const auto& nlCon : nlConstraints_) {
            // Add all variables used in nonlinear expressions
            nonlinearVars_.insert(nlCon.argVars.begin(), nlCon.argVars.end());

            // Also add result variable if it exists
            if (nlCon.resultVar >= 0) {
                nonlinearVars_.insert(nlCon.resultVar);
            }
        }

        std::vector<KNINT> linearVarIndices;
        linearVarIndices.reserve(numVars - nonlinearVars_.size());
        for (int i = 0; i < numVars; i++) {
            if (nonlinearVars_.find(i) == nonlinearVars_.end()) {
                linearVarIndices.push_back(i);
            }
        }

        // Set properties for linear variables
        if (!linearVarIndices.empty()) {
            std::vector<int> properties(linearVarIndices.size(), KN_VAR_LINEAR);

            KNITROMP_CCALL(KN_set_var_properties(
                lp(),
                linearVarIndices.size(),
                linearVarIndices.data(),
                properties.data()
            ));
            fmt::print("Marked these vars as integer: ");
            for (size_t i = 0; i < linearVarIndices.size(); ++i) {
                fmt::print("{} ", linearVarIndices[i]);
			}
			fmt::print("\n");
        }
        
        nonlinearVars_.clear();
    }


    void KnitrompModelAPI::FinishProblemModificationPhase() {
        // Register callback for nonlinear constraints
        if (nlObjective_.has_value() || !nlConstraints_.empty()) {
           
            // Note: this does not work currently, not sure why
            //markLinearVariables();
            
            // Add function evaluation callback (mandatory)
            std::vector<int> conIndices;
            for (const auto& nlCon : nlConstraints_) 
                conIndices.push_back(nlCon.knitroConIndex);
            CB_context_ptr cbContext;
            KNITROMP_CCALL(KN_add_eval_callback(
                lp(),
                nlObjective_.has_value() ? KNTRUE : KNFALSE, // if we evaluate obj
                static_cast<int>(conIndices.size()),
                conIndices.data(),
                evalNonlinearConstraint,
                &cbContext));

            // Store pointer to this object as user data
            KNITROMP_CCALL(KN_set_cb_user_params(
                lp(), cbContext, static_cast<void*>(this)));
            
            if (useJacobian)
            {
                // Add Jacobian // gradient callback
                std::vector<int> jacIndexCons;
                std::vector<int> jacIndexVars;
                std::vector<int> objGradIndexVars;
                // Build up the jacobian (and objective gradient) structure
                for (const auto& nlCon : nlConstraints_) {
                    // Add Jacobian structure for this constraint
                    for (int argVar : nlCon.argVars) {
                        jacIndexCons.push_back(nlCon.knitroConIndex);
                        jacIndexVars.push_back(argVar);
                    }
                }
                if (nlObjective_.has_value()) {
                    for (int argVar : nlObjective_->argVars) {
                        objGradIndexVars.push_back(argVar);
                    }
                }

                // Set up gradient structure for all constraints and obj
                KNITROMP_CCALL(KN_set_cb_grad(
                    lp(), cbContext,
                    objGradIndexVars.size(), objGradIndexVars.data(), // objective sparse gradient data
                    static_cast<int>(jacIndexVars.size()),
                    jacIndexCons.data(),
                    jacIndexVars.data(),
                    evalNonlinearConstraint));

            }
            
            if (useHessian) {
                // Build UNIFIED sparsity pattern to eliminate duplicates
                hessSparsityMap_.clear();

                // Helper lambda to add sparsity entries
                auto addSparsity = [this](const NonlinearConstraintData& nlData) {
                    for (const auto& [localRow, localCol] : nlData.hessianSparsity) {
                        // Convert local indices to global variable indices
                        int globalRow = nlData.argVars[localRow];
                        int globalCol = nlData.argVars[localCol];

                        // Ensure row >= col for lower triangle
                        if (globalCol > globalRow)
                            std::swap(globalRow, globalCol);

                        auto key = std::make_pair(globalRow, globalCol);

                        // Add to map if not already present
                        if (hessSparsityMap_.find(key) == hessSparsityMap_.end()) {
                            hessSparsityMap_[key] = hessSparsityMap_.size();
                        }
                    }
                    };

                // Add objective Hessian sparsity
                if (nlObjective_.has_value()) {
                    addSparsity(*nlObjective_);
                }

                // Add constraint Hessian sparsity
                for (const auto& nlCon : nlConstraints_) {
                    addSparsity(nlCon);
                }

                // Now build the mapping for each constraint/objective
                // This maps local sparsity index → global Hessian array index
                auto buildGlobalIndices = [this](NonlinearConstraintData& nlData) {
                    nlData.hessianGlobalIndices.clear();
                    nlData.hessianGlobalIndices.reserve(nlData.hessianSparsity.size());

                    for (const auto& [localRow, localCol] : nlData.hessianSparsity) {
                        int globalRow = nlData.argVars[localRow];
                        int globalCol = nlData.argVars[localCol];

                        if (globalCol > globalRow)
                            std::swap(globalRow, globalCol);

                        auto key = std::make_pair(globalRow, globalCol);
                        nlData.hessianGlobalIndices.push_back(hessSparsityMap_.at(key));
                    }
                    };

                if (nlObjective_.has_value()) {
                    buildGlobalIndices(*nlObjective_);
                }

                for (auto& nlCon : nlConstraints_) {
                    buildGlobalIndices(nlCon);
                }

                // Convert unified map to vectors for Knitro
                std::vector<int> hessIndexRows(hessSparsityMap_.size());
                std::vector<int> hessIndexCols(hessSparsityMap_.size());
                for (const auto& [coords, idx] : hessSparsityMap_) {
                    hessIndexRows[idx] = coords.first;
                    hessIndexCols[idx] = coords.second;
                }
            KN_set_cb_hess(lp(), cbContext,
                hessIndexRows.size(),
                hessIndexRows.data(),
                hessIndexCols.data(),
                evalHessian);
            } // if Use hessian

        } // if Has non linear constraints


        // Provide printing function
        set_format_model([this](fmt::MemoryWriter &w) {
            this->formatModel(w);
            });

        if (printProblem > 0) {
            fmt::MemoryWriter w;
            formatModel(w);
            fmt::print(w.str());
		}
    }


    void KnitrompModelAPI::formatModel(fmt::MemoryWriter &w) const {

        double hs;
        int direction;
        KN_get_obj_goal(lp(), &direction);

        // Vars
        int nv = NumVars();
		std::vector<double> varLobnds(nv);
		std::vector<double> varUpbnds(nv);
		std::vector<int> varTypes_(nv);
        KNITROMP_CCALL(KN_get_var_lobnds_all(lp(), varLobnds.data()));
        KNITROMP_CCALL(KN_get_var_upbnds_all(lp(), varUpbnds.data()));
		KNITROMP_CCALL(KN_get_var_types_all(lp(), varTypes_.data()));

        // Get var names
        std::vector<std::vector<char>> nameBuffers(nv, std::vector<char>(20));
        std::vector<char*> varNames(nv);
        for (int i = 0; i < nv; ++i) {
            varNames[i] = nameBuffers[i].data();
        }   
        KNITROMP_CCALL(KN_get_var_names_all(lp(), 20, varNames.data()));

        // Create lambda to resolve variable names
        auto getVarName = [&varNames](int idx, fmt::MemoryWriter &w) -> void {
            if (varNames[idx] == nullptr)
                w << "x[" << idx << "]";
            else
                w << varNames[idx];
          };
        for (int i = 0; i < nv; i++) {
            if (varNames[i] == nullptr)
                w << "x[" << i << "]";
            else
                w << varNames[i];
            w << " in [" << varLobnds[i] << ", " << varUpbnds[i] << "]";
            if (varTypes_[i] == KN_VARTYPE_INTEGER)
                w << " (integer)";
			if (varTypes_[i] == KN_VARTYPE_BINARY)
                w << " (binary)";
            w << "\n";
		}


        if (direction == KN_OBJGOAL_MAXIMIZE)
            w << "maximize ";
        else
            w << "minimize ";
        // Obj
        if (nlObjective_.has_value())
            ExpressionData::FormatExpression(w, *nlObjective_->originalExpr, getVarName);
        else
            printConstraints_.at(-1).printObjective(direction == KN_OBJGOAL_MAXIMIZE, w, getVarName);
        w << "\n";

        // Constraints
         // Get var names
        auto nc = NumCons();
        std::vector<std::vector<char>> cnameBuffers(nc, std::vector<char>(20));
        std::vector<char*> conNames(nc);
        for (int i = 0; i < nc; ++i) {
            conNames[i] = cnameBuffers[i].data();
        }
        KNITROMP_CCALL(KN_get_con_names_all(lp(), 20, conNames.data()));
        
        for (size_t i = 0; i < NumCons(); i++){
            if (conNames[i] == nullptr)
                w << "c" << i;
            else
                w << conNames[i];
            w << ": ";

            // If linear
            if (printConstraints_.find(i) != printConstraints_.end()) {
                printConstraints_.at(i).print(w, getVarName);
                continue;
			}

            // If non linear
            const auto& c = nlConstraints_[knitroToNLConstraintIndex_.at(i)];

            // LB
            KNITROMP_CCALL(KN_get_con_lobnd(lp(), c.knitroConIndex, &hs));
            if ((c.resultVar < 0) && (hs != -Infinity()))
                w << hs << "<= ";

            // or result variable
            if (c.resultVar >= 0) { 
                getVarName(c.resultVar, w);
				w << " = ";
            }
            ExpressionData::FormatExpression(w, *c.originalExpr, getVarName);
            
            // Show linear terms
            if (!c.linearIndices_.empty()) {
                for (size_t j = 0; j < c.linearIndices_.size(); j++) {
                    PrintConstraint::PrintCoef(c.linearCoeffs_[j], w, false);
                    getVarName(c.linearIndices_[j], w);
                }
            }
            // UB
            KNITROMP_CCALL(KN_get_con_upbnd(lp(), c.knitroConIndex, &hs));
            if ((c.resultVar < 0) && (hs != Infinity()))
                w << " <= " << hs;
            w << "\n";
        }

    }
    ExpressionData KnitrompModelAPI::AddExpression(const NLAffineExpression& ae) {
        ExpressionData exp;
        
        double ct = GetConstTerm(ae);
        bool hasConst = ct != 0;
        int size = GetLinSize(ae);

		for (int i = 0; i < size; ++i) {
            
            auto term = GetLinTerm(ae, i);
            auto coeff = GetLinCoef(ae, i);
            if (coeff != 1.0)
                term = ExpressionData::Multiply(ExpressionData::MakeConstantExpr(coeff), term);

            if (i == 0)
                exp = term;
            else
				exp = ExpressionData::Add(exp, term);
        }
        if (hasConst)
        {
            auto ce = ExpressionData::MakeConstantExpr(ct);
            if (size > 0)
                exp = ExpressionData::Add(exp, ce);
            else
                exp = ce;
        }
        return exp;
    }
    ExpressionData KnitrompModelAPI::AddExpression(const NLQuadExpression& qe) {
        
        ExpressionData exp;
        double ct = GetConstTerm(qe);
        bool hasConst = ct != 0;
        
        int size = GetLinSize(qe);
        for (int i = 0; i < size; ++i) {

            auto term = GetLinTerm(qe, i);
            auto coeff = GetLinCoef(qe, i);
            if (coeff != 1.0)
                term = ExpressionData::Multiply(ExpressionData::MakeConstantExpr(coeff), term);

            if (exp.isEmpty())
                exp = term;
            else
                exp = ExpressionData::Add(exp, term);
        }

        size = GetQuadSize(qe);
        for (int i = 0; i < size; ++i) {
            auto x1 = GetQuadTerm1(qe, i);
			auto x2 = GetQuadTerm2(qe, i);
            auto term = ExpressionData::Multiply(x1, x2);
            auto coeff = GetQuadCoef(qe, i);
            if (coeff != 1.0)
                term = ExpressionData::Multiply(ExpressionData::MakeConstantExpr(coeff), term);
            if (exp.isEmpty())
                exp = term;
            else
                exp = ExpressionData::Add(exp, term);
        }

        // Constant term
        if (ct && (exp.isEmpty()))
            exp = ExpressionData::MakeConstantExpr(ct);
        else
            if(ct != 0.0)
			    exp = ExpressionData::Add(exp, ExpressionData::MakeConstantExpr(ct));
        return exp;


    }

    
    // var >= expr.
    void KnitrompModelAPI::AddConstraint(const NLAssignGE& nl) {
        auto resultVar = GetVariable(nl);
        auto exp = GetExpression(nl);
        
        // AddAssign... will create exp - resultVar <SENSE> 0
		// To create var >= expr, we have to rewrite it as expr - var <= 0
		addAssignNonlinearConstraint(LE, 0.0, resultVar, exp, nl.name());
    }
    // var <= expr.
    void KnitrompModelAPI::AddConstraint(const NLAssignLE& nl) {
        auto resultVar = GetVariable(nl);
        auto exp = GetExpression(nl);
        // AddAssign... will create exp - resultVar <SENSE> 0
        // To create var <= expr, we have to rewrite it as expr - var >= 0
        addAssignNonlinearConstraint(GE, 0.0, resultVar, exp, nl.name());
    }

    void KnitrompModelAPI::AddConstraint(const NLAssignEQ& nl) {
        auto resultVar = GetVariable(nl);
        auto exp = GetExpression(nl);
        addAssignNonlinearConstraint(EQ, 0.0, resultVar, exp, nl.name());
	}

    void KnitrompModelAPI::SetNLObjective(int i, const NLObjective& nl) {
        auto exp = GetExpression(nl);

        if (nl.obj_sense() == obj::MAX) {
            KNITROMP_CCALL(KN_set_obj_goal(lp(), KN_OBJGOAL_MAXIMIZE));
        }
        const auto &lin = nl.GetLinTerms();
        if (lin.size() > 0) {
            exp.reserveLinear(lin.size());
            for (int i = 0; i < lin.size(); ++i)
				exp.addLinear(lin.pvars()[i], lin.pcoefs()[i]);

        }

        // Add linear part
        if (exp.linearCoeffs().size() > 0)
            KNITROMP_CCALL(KN_add_obj_linear_struct(lp(), exp.linearCoeffs().size(), 
                exp.linearIndices().data(), exp.linearCoeffs().data()));

        storeNonLinearData(-1, exp, -1, nl.name());
        
    }
    void KnitrompModelAPI::AddConstraint(const NLConstraint& nl) {
        auto exp = GetExpression(nl);
        double lhs = GetLower(nl), rhs = GetUpper(nl);
        char type;
        double range;
        double* prange = nullptr;

        if (GetLinSize(nl) > 0) {
            exp.reserveLinear(GetLinSize(nl));
            for (int i = 0; i < GetLinSize(nl); ++i)
                exp.addLinear(GetLinVar(nl, i), GetLinCoef(nl, i));
        }
        
        // Add a nonlinear constraint to Knitro: resultVar - f(argVars) = 0
        int conIndex;
        KNITROMP_CCALL(KN_add_con(lp(), &conIndex));

        lhs = lhs < -Infinity() ? -Infinity() : lhs;
        rhs = rhs > Infinity() ? Infinity() : rhs;

        if (lhs == rhs)
            KNITROMP_CCALL(KN_set_con_eqbnd(lp(), conIndex, lhs));
        else {
            if (lhs > -Infinity())
                KNITROMP_CCALL(KN_set_con_lobnd(lp(), conIndex, lhs));
            if (rhs < Infinity())
                KNITROMP_CCALL(KN_set_con_upbnd(lp(), conIndex, rhs));
        }
        // Add linear part
        if (exp.linearCoeffs().size() > 0)
            KNITROMP_CCALL(KN_add_con_linear_struct_one(lp(), exp.linearCoeffs().size(), conIndex,
                exp.linearIndices().data(), exp.linearCoeffs().data()));

        storeNonLinearData(conIndex, exp, -1, nl.name());

    }

    template <class MPExpr> ExpressionData KnitrompModelAPI::AddUnaryExpression(ExpressionData::OpType type,
        const MPExpr& expr) {
        return ExpressionData::MakeUnaryExpr(type,
            GetArgExpression(expr, 0));
    }


    ExpressionData KnitrompModelAPI::AddExpression(const SinExpression& ex) {
        return AddUnaryExpression(ExpressionData::SIN, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const CosExpression& ex) {
        return AddUnaryExpression(ExpressionData::COS, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const TanExpression& ex) {
        return AddUnaryExpression(ExpressionData::TAN, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const SinhExpression& ex) {
        return AddUnaryExpression(ExpressionData::SINH, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const CoshExpression& ex) {
        return AddUnaryExpression(ExpressionData::COSH, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const TanhExpression& ex) {
        return AddUnaryExpression(ExpressionData::TANH, ex);
    }


    ExpressionData KnitrompModelAPI::AddExpression(const AsinExpression& ex) {
        return AddUnaryExpression(ExpressionData::ASIN, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const AcosExpression& ex) {
        return AddUnaryExpression(ExpressionData::ACOS, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const AtanExpression& ex) {
        return AddUnaryExpression(ExpressionData::ATAN, ex);
    }

    ExpressionData KnitrompModelAPI::AddExpression(const AsinhExpression & ex) {
        return AddUnaryExpression(ExpressionData::ASINH, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const AcoshExpression & ex) {
        return AddUnaryExpression(ExpressionData::ACOSH, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const AtanhExpression& ex) {
        return AddUnaryExpression(ExpressionData::ATANH, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const LogExpression& ex) {
        return AddUnaryExpression(ExpressionData::LOG, ex);
    }
    ExpressionData KnitrompModelAPI::AddExpression(const LogAExpression& ex) {
        auto log_x = ExpressionData::MakeUnaryExpr(ExpressionData::LOG, GetArgExpression(ex, 0));
        double base = GetParameter(ex, 0);
        double log_base = std::log(base); 
        auto log_base_expr = ExpressionData::MakeConstantExpr(log_base);
        return ExpressionData::BinaryOp(ExpressionData::DIV, log_x, log_base_expr);
    }

    ExpressionData KnitrompModelAPI::AddExpression(const DivExpression& e) {
        ExpressionData exp;
        return exp.BinaryOp(ExpressionData::DIV, GetArgExpression(e, 0),
            GetArgExpression(e, 1));
    }

    ExpressionData KnitrompModelAPI::AddExpression(const PowExpression& e) {
        ExpressionData exp;
        return exp.BinaryOp(ExpressionData::POW, GetArgExpression(e, 0),
            GetArgExpression(e, 1));
    }

    ExpressionData KnitrompModelAPI::AddExpression(const PowConstExpExpression& e) {
        ExpressionData exp;
        return exp.BinaryOp(ExpressionData::POW, GetArgExpression(e, 0),
            ExpressionData::MakeConstantExpr(GetParameter(e, 0)));
    }

    CppAD::AD<double> buildAD(const ExpressionData& d, const std::map<int, CppAD::AD<double>>& varMap) {
        switch (d.op()) {
        case ExpressionData::VAR:
            return varMap.at(d.varIndex());
        case ExpressionData::CONST:
            return CppAD::AD<double>(d.constValue());
        case ExpressionData::ADD:
            return buildAD(*d.left(), varMap) + buildAD(*d.right(), varMap);
        case ExpressionData::SUB:
            return buildAD(*d.left(), varMap) - buildAD(*d.right(), varMap);
        case ExpressionData::MUL:
            return buildAD(*d.left(), varMap) * buildAD(*d.right(), varMap);
        case ExpressionData::DIV:
            return buildAD(*d.left(), varMap) / buildAD(*d.right(), varMap);
        case ExpressionData::SIN:
            return CppAD::sin(buildAD(*d.left(), varMap));
        case ExpressionData::COS:
            return CppAD::cos(buildAD(*d.left(), varMap));
        case ExpressionData::TAN:
            return CppAD::tan(buildAD(*d.left(), varMap));
        case ExpressionData::LOG:
            return CppAD::log(buildAD(*d.left(), varMap));
        case ExpressionData::EXP:
            return CppAD::exp(buildAD(*d.left(), varMap));
        case ExpressionData::POW:
            return CppAD::pow(buildAD(*d.left(), varMap), buildAD(*d.right(), varMap));
        case ExpressionData::SQRT:
            return CppAD::sqrt(buildAD(*d.left(), varMap));
        case ExpressionData::ASIN:
            return CppAD::asin(buildAD(*d.left(), varMap));
        case ExpressionData::ACOS:
            return CppAD::acos(buildAD(*d.left(), varMap));
        case ExpressionData::ATAN:
            return CppAD::atan(buildAD(*d.left(), varMap));
        case ExpressionData::SINH:
            return CppAD::sinh(buildAD(*d.left(), varMap));
        case ExpressionData::COSH:
            return CppAD::cosh(buildAD(*d.left(), varMap));
        case ExpressionData::TANH:
            return CppAD::tanh(buildAD(*d.left(), varMap));
        case ExpressionData::ASINH:
            return CppAD::asinh(buildAD(*d.left(), varMap));
        case ExpressionData::ACOSH:
            return CppAD::acosh(buildAD(*d.left(), varMap));
        case ExpressionData::ATANH:
            return CppAD::atanh(buildAD(*d.left(), varMap));
        default:
            throw std::runtime_error("Unknown operation type");
        }
    }

    // Create a CppAD tape from this expression tree
    std::pair<CppAD::ADFun<double>, std::vector<int>> createTape(const ExpressionData& d) {
        // Extract variable indices in sorted order

        std::vector<int> varIndices(d.usedVars().begin(), d.usedVars().end());

        // Create independent variables
        std::vector<CppAD::AD<double>> X(varIndices.size());
        CppAD::Independent(X);

        // Build variable mapping
        std::map<int, CppAD::AD<double>> varMap;
        for (size_t i = 0; i < varIndices.size(); i++) {
            varMap[varIndices[i]] = X[i];
        }

        // Build the CppAD expression from AST
        CppAD::AD<double> result = buildAD(d, varMap);

        // Create tape
        std::vector<CppAD::AD<double>> Y(1);
        Y[0] = result;
        CppAD::ADFun<double> tape(X, Y);

        return { std::move(tape), std::move(varIndices) };
    }

} // namespace mp
