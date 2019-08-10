#ifndef _MATH_LINK_TABLE_EVALUATOR_BASE_2009_10_17_H_
  #define _MATH_LINK_TABLE_EVALUATOR_BASE_2009_10_17_H_

  #include <string>
  #include <vector>

  class MathLinkBase;

  class MathLinkTableEvaluatorBase
  {
  private:

    const MathLinkTableEvaluatorBase& operator=(const MathLinkTableEvaluatorBase&);
    MathLinkTableEvaluatorBase(const MathLinkTableEvaluatorBase&);

  protected:

    const MathLinkBase& my_ML;
    const std::string   my_expression;

    mutable std::vector<std::vector<std::string> > str_values;

  protected:

    MathLinkTableEvaluatorBase(const MathLinkBase& ml,
                               const std::string&  expr) : my_ML        (ml),
                                                           my_expression(expr) { }

    static void Replace_StarHat_With_e(std::string& str);

  public:

    virtual ~MathLinkTableEvaluatorBase() { }

    const std::vector<std::vector<std::string> >& StringValues(void) const { return str_values; }

    virtual bool Evaluate(void) const = 0;
  };
  
#endif // _MATH_LINK_TABLE_EVALUATOR_BASE_2009_10_17_H_
