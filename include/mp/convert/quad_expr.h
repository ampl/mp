#ifndef QUAD_EXPR_H
#define QUAD_EXPR_H

#include "mp/convert/affine_expr.h"

namespace mp {

class QuadTerms {
private:
 class Term {
  private:
   using VarPair = std::pair<int, int>;
   VarPair vars_;
   double coef_;

   friend class QuadTerms;

   Term(VarPair vars, double coef) : vars_(vars), coef_(coef) {}

  public:
   VarPair vars() const { return vars_; }
   double coef() const { return coef_; }
   void set_coef(double c) { coef_=c; }
 };
 std::vector<Term> terms_;

public:
 QuadTerms() { }

 bool empty() const { return terms_.empty(); }
 int num_terms() const { return static_cast<int>(terms_.size()); }
 int capacity() const { return static_cast<int>(terms_.capacity()); }

 Term::VarPair vars(int i) const { return terms_[i].vars(); }
 double coef(int i) const { return terms_[i].coef(); }
 void set_coef(int i, double c) { terms_[i].set_coef(c); }

 typedef std::vector<Term>::const_iterator const_iterator;

 const_iterator begin() const { return terms_.begin(); }
 const_iterator end() const { return terms_.end(); }

 typedef std::vector<Term>::iterator iterator;

 iterator begin() { return terms_.begin(); }
 iterator end() { return terms_.end(); }

 void AddTerm(Term::VarPair vars, double coef) {
   terms_.push_back(Term(vars, coef));
 }

 void AddTerms(const QuadTerms& li) {
   terms_.insert(end(), li.begin(), li.end());
 }

 void Reserve(int num_terms) {
   terms_.reserve(num_terms);
 }

 /// Arithmetic
 void Negate() {
   for (auto& term: *this)
     term.set_coef(-term.coef());
 }

 void Add(const QuadTerms& ae) {
   this->Reserve(this->num_terms() + ae.num_terms());
   this->AddTerms(ae); // eliminate duplicates when?
 }

 void Subtract(QuadTerms&& ae) {
   ae.Negate();
   Add(ae);
 }

};

class QuadExpr {
public:
  QuadExpr() { }
  QuadExpr(AffineExpr&& ae) : ae_(std::move(ae)) { }

  using Constant = AffineExpr::Constant;
  using Variable = AffineExpr::Variable;

  /// Getters
  bool is_constant() const { return is_affine() && GetAE().is_constant(); }
  double constant_term() const { return GetAE().constant_term(); }
  bool is_variable() const { return is_affine() && GetAE().is_variable(); }
  int get_representing_variable() const {
    assert(is_affine());
    return GetAE().get_representing_variable();
  }

  bool is_affine() const { return GetQT().empty(); }
  const AffineExpr& GetAE() const { return ae_; }
  AffineExpr& GetAE() { return ae_; }
  bool is_quadratic() const { return !is_affine(); }
  const QuadTerms& GetQT() const { return qt_; }
  QuadTerms& GetQT() { return qt_; }


  /// Modifiers
  void constant_term(double v) { GetAE().constant_term(v); }
  void add_to_constant(double a) { GetAE().add_to_constant(a); }
  void AddLinearTerm(int var_index, double coef) {
    GetAE().AddTerm(var_index, coef);
  }

  void Add(const QuadExpr& qe) {
    GetAE().Add(qe.GetAE());
    GetQT().Add(qe.GetQT());
  }

  void Subtract(QuadExpr&& qe) {
    GetAE().Subtract(std::move(qe.GetAE()));
    GetQT().Subtract(std::move(qe.GetQT()));
  }

  void Negate() {
    GetAE().Negate();
    GetQT().Negate();
  }
private:
  AffineExpr ae_;
  QuadTerms qt_;
};

}

#endif // QUAD_EXPR_H
