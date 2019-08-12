1.  b_0 is a card-carrying PITA. I invariably have to subtract it off.
    Change this mini-lib - quickly - to start with index 1 in the
    numerator function.

2.  Actually apply the tail function!
    This involves stopping the loop at an even iteration.
    Then do one "loop" with w replacing something.
    Evaluate w and if nonzero...

3.  [Done] Allow setting of number of iterations and tolerance in
          all the classes.

          I just moved some data members to public and made them non-const.

4.  Try to make constexpr.  Think of different way to handle errors.

5.  Finish forward cfrac!

6.  Make an exception type that helps people try something different:
    max_iters, epsilon, result, previous delta.

7.  Learn how to LaTex K_{m=1}^{\infty}.
