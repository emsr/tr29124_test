#ifndef _MATH_LINK_TABLE_EVALUATOR_COMPLEX_2009_12_18_H_
  #define _MATH_LINK_TABLE_EVALUATOR_COMPLEX_2009_12_18_H_

  #include <TableEvaluator/MathLinkTableEvaluatorBase.h>

  class MathLinkTableEvaluatorComplex : public MathLinkTableEvaluatorBase
  {
  public:

    MathLinkTableEvaluatorComplex(const MathLinkBase& ml,
                                  const std::string&  expr) : MathLinkTableEvaluatorBase(ml, expr) { }

    virtual ~MathLinkTableEvaluatorComplex() { }

  private:

    virtual bool Evaluate(void) const;
  };
  
#endif // _MATH_LINK_TABLE_EVALUATOR_COMPLEX_2009_12_18_H_
