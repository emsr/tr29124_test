#ifndef _MATH_LINK_TABLE_EVALUATOR_2009_10_17_H_
  #define _MATH_LINK_TABLE_EVALUATOR_2009_10_17_H_

  #include <TableEvaluator/MathLinkTableEvaluatorBase.h>

  class MathLinkTableEvaluator : public MathLinkTableEvaluatorBase
  {
  public:

    MathLinkTableEvaluator(const MathLinkBase& ml,
                           const std::string&  expr) : MathLinkTableEvaluatorBase(ml, expr) { }

    virtual ~MathLinkTableEvaluator() { }

  private:

    virtual bool Evaluate(void) const;
  };
  
#endif // _MATH_LINK_TABLE_EVALUATOR_2009_10_17_H_
