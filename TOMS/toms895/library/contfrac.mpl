##### ========================================
##### ===== CONTINUED FRACTION SUBMODULE =====
##### ========================================

contfrac := module()

  description "CFSF continued fraction submodule";

  export
    argtypes, argnames, reqargs,
    create, formula,
    nthpelement, nthpnumer, nthpdenom,
    nthnumer, nthdenom, nthapprox,
    nthtail, tailestimate,
    transform;

  local
    setup, cleanup;

  option
    load = setup, unload = cleanup;

##### ----- LOADING AND UNLOADING -----

  setup := proc()

    description "setup procedure for CFSF[contfrac]";

    global `type/simregular`;

    USERINFO( 5, ':-CFSF', `setup procedure for CFSF[contfrac]` );

    `type/simregular` := proc( aformula )

      local tformula, t;

      tformula := eval( aformula );

      if type( tformula, table ) and validate( tformula ) and tformula['type'] = 'CFSF:-contfrac' then

        t := { seq( tformula['general'][j][2], j=1..nops(tformula['general']) ) };

        if assigned( tformula['begin'] ) then
          t := t union { seq( tformula['begin'][i][2], i=1..nops(tformula['begin']) ) };
        end if;

        if nops( t ) = 1 and t[1] = 1 then
          return true;
        end if;

      end if;

      return false;

    end proc;

    NULL

  end proc;

  cleanup := proc()

    description "cleanup procedure for CFSF[contfrac]";

    USERINFO( 5, ':-CFSF', `cleanup procedure for CFSF[contfrac]` );

    # `type/simregular` := evaln( `type/simregular` );

    NULL

  end proc;

##### ----- CREATING AND VALIDATING FORMULAS -----

  argtypes := table([ op(op(op(cargtypes))),
    'begin' = list([exprtype,exprtype]),
    'even' = [exprtype,exprtype],
    'odd' = [exprtype,exprtype],
    'general' = list([exprtype,exprtype]),
    'modification' = exprtype
  ]);

  argnames := map( op, { indices(argtypes) } );

  reqargs := { op(reqcargs), 'general' };

  create := proc( newformula::table )

    if not assigned( newformula['index'] ) then
      newformula['index'] := 'm';
    end if;

  end proc;

##### ----- RETRIEVING FORMULAS -----

  formula := proc( tformula::table )

    local tmodification;

    # if a modification argument was specified, we add it to formula

    if nargs = 2 then
      if type( args[2], equation ) and op(1,args[2]) = 'modification' then
        tmodification := op(2,args[2]);
      elif not type( args[2], set ) then
        ERROR( "`formula' expects its 2nd argument to be of type `set' or to be of the form `modification = <expression>', but got `%1'", args[2] );
      end if;
    elif nargs = 3 then
      if not type( args[2], set ) then
        ERROR( "`formula' expects its 2nd argument to be of type `set', but got `%1'", args[2] );
      elif type( args[3], equation ) and op(1,args[3]) = 'modification' then
        tmodification := op(2,args[3]);
      else
        ERROR( "`formula' expects its 3rd argument to be of the form `modification = <expression>', but got `%1'", args[3] );
      end if;
    end if;

    if assigned( tmodification ) then
      if assigned( tformula['modification'] ) then
        WARNING( "the formula previously had value `%1' for argument `modification'", eval(tformula['modification']) );
        WARNING( "replacing with new value `%1'", eval(tmodification) );
      end if;
      tformula['modification'] := tmodification;
    end if;

  end proc;

##### ----- RETRIEVING PARTS OF FORMULAS -----

  nthpelement := proc( tformula::table, n::nonnegint )

    # NOTE : we use an extra eval so that assumptions show up in the output

    if n = 0 then

      if assigned( tformula['front'] ) then
        eval( tformula['front'] );
      else
        0;
      end if;

    elif assigned( tformula['begin'] ) and n <= nops( tformula['begin'] ) then

      if n = 1 and assigned( tformula['factor'] ) then
        eval( [ tformula['factor'] * tformula['begin'][1][1], tformula['begin'][1][2] ] );
      else
        eval( tformula['begin'][n] );
      end if;

    else

      if n = 1 and assigned( tformula['factor'] ) then
        eval( [ tformula['factor'] * subs( tformula['index']=1, tformula['general'][1][1] ), subs( tformula['index']=1, tformula['general'][1][2] ) ] );
      elif not assigned( tformula['begin'] ) then
        eval( subs( tformula['index']=n, tformula['general'][(n-1) mod nops(tformula['general'])+1] ) );
      else
        eval( subs( tformula['index']=n, tformula['general'][(n-nops(tformula['begin'])-1) mod nops(tformula['general'])+1] ) );
      end if;

    end if;

  end proc;

  nthpnumer := proc( tformula::table, n::nonnegint )

    nthpelement(tformula,n)[1];

  end proc;

  nthpdenom := proc( tformula::table, n::nonnegint )

    if n = 0 then
      1;
    else
      nthpelement(tformula,n)[2];
    end if;

  end proc;

##### ----- EVALUATING FORMULAS -----

  nthapprox := proc( tformula::table, n::nonnegint, tsubs::table )

    local tsubstitutions, tmodification, expr, ai, bi, i;

    tsubstitutions, tmodification := constructsubs( tformula, tsubs );

    # if assigned, we start by inserting the modification
    # and set F_{n+1} = w, otherwise F_{n+1} = 0

    # NOTE : since tmodification can be the result of a call to
    # tailestimate (in constructsubs), we cannot use subs here because
    # it also replaces tformula['index'] inside formula objects,
    # making the evaluation invalid; algsubs does not replace inside
    # function calls

    if not tmodification = 0 then
      # USERINFO( 4, ':-CFSF', `modification argument equals`, tmodification );
      expr := subs( tsubstitutions, eval( algsubs( tformula['index'] = n, tmodification ) ) );
    else
      expr := 0;
    end if;

    # using backward recurrence : F_{i} = a_{i} / ( b_{i} + F_{i+1} )
    # for i = n,n-1,..,1 and finally F_{0} = b_{0} + F_{1}

    # NOTE : we use try-catch block here to avoid numeric exceptions when
    # executing the backward recurrence (e.g., division by zero, which
    # occurred in nthapprox(formula("EF.exp.cfrac.01"),5,z=15));
    # if such an exception does occur, we first do a symbolic evaluation
    # and substitute the values afterwards

    # NOTE : we normalize after every 10 steps, since otherwise the
    # final normalization needs too much memory (and time as well)

    try

      for i from n by -1 to 1 do
        ai, bi := op( subs( tsubstitutions, nthpelement( tformula, i ) ) );
        expr := ai / ( bi + expr );
        if assigned( tsubs['variable'] ) and type( tsubs['variable'], complex( float ) ) then
          expr := evalf( expr );
        elif defaults[evaluate_simplify] > 0 and i mod defaults[evaluate_simplify] = 0 then
          expr := normal( expr );
        end if;
      end do;

      expr := subs( tsubstitutions, nthpelement( tformula, 0 ) ) + expr;

    catch "numeric exception" :

      if assigned( tsubs['variable'] ) then
        WARNING( "encountered `%1', using symbolic method instead", lasterror );
        expr := subs( tsubstitutions, :-CFSF:-nthapprox( tformula, n ) );
      else
        # WARNING( "encountered `%1'", lasterror );
        ERROR( lasterror );
      end if;

    end try;

    # if the variable is of type float, then we do an extra evalf
    # to eliminate constants from the output (e.g. Pi)

    if assigned( tsubs['variable'] ) and type( tsubs['variable'], complex(float) ) then
      expr := evalf( expr );
    end if;

    # finally, we normalize and sort the output
    # (since this is not done automatically)

    return sort( normal( expr ) );

  end proc;

  nthnumer := proc( tformula::table, n::nonnegint, tsubs::table )

    description "function to compute the nth numerator of a continued fraction";

    local tsubstitutions, tmodification, Aim2, Aim1, Ai, ai, bi, i;

    tsubstitutions, tmodification := constructsubs( tformula, tsubs );

    # using recurrence relation : A_{i} = b_{i} * A_{i-1} + a_{i} * A_{i-2}
    # for i = 1,2,..,n with initial conditions A_{-1} = 1 and A_{0} = b_{0}

    # NOTE : we expand at every step to avoid problems
    # with the final normalization

    Aim1 := 1;
    Ai := subs( tsubstitutions, nthpelement( tformula, 0 ) );

    for i from 1 to n do
      Aim2 := Aim1;
      Aim1 := Ai;
      ai, bi := op( subs( tsubstitutions, nthpelement( tformula, i ) ) );
      Ai := expand( bi * Aim1 + ai * Aim2 );
      if assigned( tsubs['variable'] ) and type( tsubs['variable'], complex( float ) ) then
        Ai := evalf( Ai );
      elif defaults[evaluate_simplify] > 0 and i mod defaults[evaluate_simplify] = 0 then
        Ai := normal( Ai );
      end if;
    end do;

    # if assigned, the modified numerator is given by A_{n} + A_{n-1} * w

    if not tmodification = 0 then
      Ai := Ai + Aim1 * subs( tsubstitutions union { tformula['index'] = n }, tmodification );
    end if;

    if assigned( tsubs['variable'] ) and type( tsubs['variable'], complex(float) ) then
      Ai := evalf( Ai );
    end if;

    return sort( normal( Ai ) );

  end proc;

  nthdenom := proc( tformula::table, n::nonnegint, tsubs::table )

    description "function to compute the nth denominator of a continued fraction";

    local tsubstitutions, tmodification, Bim2, Bim1, Bi, ai, bi, i;

    tsubstitutions, tmodification := constructsubs( tformula, tsubs );

    # using recurrence relation : B_{i} = b_{i} * B_{i-1} + a_{i} * B_{i-2}
    # for i = 1,2,..,n with initial conditions B_{-1} = 0 and B_{0} = 1

    Bim1 := 0;
    Bi := 1;

    for i from 1 to n do
      Bim2 := Bim1;
      Bim1 := Bi;
      ai, bi := op( subs( tsubstitutions, nthpelement( tformula, i ) ) );
      Bi := expand( bi * Bim1 + ai * Bim2 );
      if assigned( tsubs['variable'] ) and type( tsubs['variable'], complex( float ) ) then
        Bi := evalf( Bi );
      elif defaults[evaluate_simplify] > 0 and i mod defaults[evaluate_simplify] = 0 then
        Bi := normal( Bi );
      end if;
    end do;

    # if assigned, the modified denominator is given by B_{n} + B{n-1} * w

    if not tmodification = 0 then
      Bi := Bi + Bim1 * subs( tsubstitutions union { tformula['index'] = n }, tmodification );
    end if;

    if assigned( tsubs['variable'] ) and type( tsubs['variable'], complex(float) ) then
      return evalf( Bi );
    end if;

    return sort( normal( Bi ) );

  end proc;

##### ----- WORKING WITH TAILS -----

  nthtail := proc( tformula::table, n::nonnegint )

    local newformula, t;

    # first we copy some arguments from the original formula; note that
    # not all arguments are copied, for instance 'labels' and 'function'
    # are not copied, since they no longer hold for the newly created
    # formula

    newformula := table();
    newformula['type'] := tformula['type'];
    newformula['index'] := tformula['index'];
    newformula['variable'] := tformula['variable'];

    if assigned( tformula['parameters'] ) then
      newformula['parameters'] := tformula['parameters'];
    end if;

    # 'front' never reappears and 'factor' only survives when n = 0

    if n = 0 and assigned( tformula['factor'] ) then
      newformula['factor'] := tformula['factor'];
    end if;

    # now adjust the 'begin' and 'general' arguments
    # depending on the value of n

    t := `if`( assigned( tformula['begin'] ), nops( tformula['begin'] ), 0 );
    if n < t then
      newformula['begin'] := tformula['begin'][n+1..-1];
      newformula['general'] := tformula['general'];
    else
      t := (n-t) mod nops( tformula['general'] ) + 1;
      if t = 1 then
        newformula['general'] := tformula['general'];
      else
        newformula['general'] := [ op(tformula['general'])[t..-1], op(tformula['general'])[1..t-1] ];
      end if;
    end if;

    newformula['general'] := subs( newformula['index'] = newformula['index'] + n, newformula['general'] );

    # if assigned, we copy the modification of the original formula

    if assigned( tformula['modification'] ) then
      newformula['modification'] := tformula['modification'];
    end if;

    # we do a safety check to make sure the formula validates correctly;
    # if not, we have an implementation problem (call it a bug)

    if not validate( newformula ) then
      ERROR( "validation failure, unable to create the tail successfully" );
    end if;

    # EXTRA : if the original formula has a label, then we add a comment
    # stating how the newly created formula was constructed from the
    # original one

    if assigned( tformula['label'] ) then
      newformula['comment'] := cat( convert( n, string ), "-th tail of continued fraction formula with label \"", tformula['label'], "\"" );
    end if;

    # here too we wrap the table in a 'formula' expression
    # so that the `print/formula' procedure will be called
    # to format the output

    return ':-CFSF:-formula'( copy(newformula) );

  end proc;

  tailestimate := proc( aformula::uneval )

    local tformula, t, tflimits, a, w, rflimits, r;

    tformula := eval( transform( aformula, 'simregular' ) );
    t := map2( op, 1, tformula['general'] );
    tflimits := convert( map( limit, t, tformula['index'] = infinity ), set );

    USERINFO( 4, ':-CFSF', `the formula has the limits`, tflimits );

    # check whether we can automatically compute a modification
    # (cfr. section 7.7 from the book); for this, the limit has
    # to be unique for all general elements

    if nops(tflimits) = 1 then

      # NOTE : when called from constructsubs, we include the substitutions
      # (passed in the last argument) when evaluating the limit

      a := eval( tflimits[1], `if`( type( args[-1], set ), args[-1], {} ) );

      if type( a, Or( infinity, SymbolicInfinity ) ) assuming additionally, op( aformula[constraints] ) then

        USERINFO( 4, ':-CFSF', `the unique limit is of type pos_infinity` );

        if type( aformula, 'simregular' ) then
          w := ( sqrt( 4 * ':-CFSF:-nthpnumer'( aformula, tformula['index'] + 1 ) + 1 ) - 1 ) / 2;
        else
          w := ( sqrt( 4 * ':-CFSF:-nthpnumer'( aformula, tformula['index'] + 1 ) / ( ':-CFSF:-nthpdenom'( aformula, tformula['index'] ) * ':-CFSF:-nthpdenom'( aformula, tformula['index'] + 1 ) ) + 1 ) - 1 ) / 2;
        end if;

      elif not coulditbe( a in RealRange( -infinity, Open(-1/4) ) ) assuming additionally, op( aformula[constraints] ) then

        USERINFO( 4, ':-CFSF', `the unique limit is of type finite` );

        a := tflimits[1];  # use the symbolic limit again
        w := ( sqrt( 4 * a + 1 ) - 1 ) / 2;

        # if requested, try to improve the previously computed modification w

        if nargs >= 2 and args[2] = 'improved' then

          rflimits := [ seq( (t[j]-a) / (t[j-1]-a), j=2..nops(t) ), subs( tformula['index'] = tformula['index'] + nops(t), t[1]-a ) / (t[-1]-a) ];
          rflimits := convert( map( limit, rflimits, tformula['index'] = infinity ), set );

          USERINFO( 4, ':-CFSF', `the formula has the improved limits`, rflimits );

          if nops(rflimits) = 1 then # and type( r, complex(numeric) ) then

            r := eval( rflimits[1] );

            USERINFO( 4, ':-CFSF', `using the improved modification as requested` );

            if type( aformula, 'simregular' ) then
              w := w + ( ':-CFSF:-nthpnumer'( aformula, tformula['index'] + 1 ) - a ) / ( 1 + ( r + 1 ) * w );
            else
              w := w + ( ':-CFSF:-nthpnumer'( subs( tformula['index'] = _||(eval(tformula['index'])), transform( aformula, 'simregular' ) ), tformula['index'] + 1 ) - a ) / ( 1 + ( r + 1 ) * w );
            end if;

          end if;

        end if;

      else
        USERINFO( 4, ':-CFSF', `cannot compute a modification for this unique limit` );
        WARNING( "a modification could not be computed, returning 0 instead" );
        return 0;
      end if;

    else
      USERINFO( 4, ':-CFSF', `the limit is not unique` );
      WARNING( "a modification could not be computed, returning 0 instead" );
      return 0;
    end if;

    if not type( aformula, 'simregular' ) then
      w := ':-CFSF:-nthpdenom'( aformula, tformula['index'] ) * w;
    end if;

    return w;

  end proc;

##### ----- TRANSFORMATIONS -----

  transform := proc( aformula::uneval, form::name )

    local tformula, newformula, t, nopst, a, b, i, j, k;

    tformula := eval( aformula );

    # first we copy some arguments from the original formula

    newformula := table();
    newformula['type'] := tformula['type'];
    newformula['index'] := tformula['index'];
    newformula['variable'] := tformula['variable'];

    if assigned( tformula['parameters'] ) then
      newformula['parameters'] := tformula['parameters'];
    end if;

    if assigned( tformula['constraints'] ) then
      newformula['constraints'] := tformula['constraints'];
    end if;

    if assigned( tformula['lhs'] ) then
      newformula['lhs'] := tformula['lhs'];
    end if;

    if assigned( tformula['function'] ) then
      newformula['function'] := tformula['function'];
    end if;

    ##### ----- even contraction -----

    if form = 'even_contraction' then

      # first we construct the front and begin parts

      t := [ seq( nthpelement( tformula, i ), i=0..nops(tformula['begin'])+2 ) ];
      a := [ seq( t[j][1], j=2..nops(t) ) ];
      b := [ seq( t[j][2], j=2..nops(t) ) ];

      if assigned( tformula['front'] ) then
        newformula['front'] := tformula['front'];
      end if;

      newformula['begin'] := [ [ a[1] * b[2], a[2] + b[1] * b[2] ] ];
      if assigned( tformula['begin'] ) then
        i := 2;
        while ( 2*i-2 <= nops(tformula['begin']) ) do
          newformula['begin'] := [ op(newformula['begin']), [ - a[2*i-2] * a[2*i-1] * b[2*i] / b[2*i-2], a[2*i] + b[2*i-1] * b[2*i] + a[2*i-1] * b[2*i] / b[2*i-2] ] ];
          i := i + 1;
        end do;
      end if;

      # then we create the part with general elements

      t := tformula['general'];  nopst := nops(t);
      a := [ seq( t[j][1], j=1..nopst ) ];
      b := [ seq( t[j][2], j=1..nopst ) ];

      j := 0;
      i := newformula['index'];
      newformula['general'] := [];

      if not assigned( tformula['begin'] ) or type( nops(tformula['begin']), even ) then
        j := j + 1;
      end if;

      for k from 1 by 1 to `if`( type( nopst, even ), nopst/2, nopst ) do
        newformula['general'] := [ op(newformula['general']), [ - subs(i=2*i-2,a[j mod nopst + 1]) * subs(i=2*i-1,a[(j+1) mod nopst + 1]) * subs(i=2*i,b[(j+2) mod nopst + 1]) / subs(i=2*i-2,b[j mod nopst + 1]), subs(i=2*i,a[(j+2) mod nopst + 1]) + subs(i=2*i-1,b[(j+1) mod nopst + 1]) * subs(i=2*i,b[(j+2) mod nopst + 1]) + subs(i=2*i-1,a[(j+1) mod nopst + 1]) * subs(i=2*i,b[(j+2) mod nopst + 1]) / subs(i=2*i-2,b[j mod nopst + 1]) ] ];
        j := (j+1) mod nopst + 1;
      end do;

      if assigned( tformula['label'] ) then
        newformula['comment'] := cat( "even contraction of continued fraction formula with label \"", tformula['label'], "\"" );
    end if;

    ##### ----- odd contraction -----

    elif form = 'odd_contraction' then

      # first we construct the front and begin parts

      t := [ seq( nthpelement( tformula, i ), i=0..nops(tformula['begin'])+5 ) ];
      a := [ seq( t[j][1], j=2..nops(t) ) ];
      b := [ seq( t[j][2], j=2..nops(t) ) ];

      if assigned( tformula['front'] ) then
        newformula['front'] := ( a[1] + tformula['front'] * b[1] ) / b[1];
      else
        newformula['front'] := a[1] / b[1];
      end if;

      newformula['begin'] := [ [ - a[1] * a[2] * b[3] / b[1], b[1] * ( a[3] + b[2] * b[3] ) + a[2] * b[3] ] ];

      if assigned( tformula['begin'] ) then
        i := 2;
        while ( 2*i-3 <= nops(tformula['begin']) ) do
          newformula['begin'] := [ op(newformula['begin']), [ - a[2*i-1] * a[2*i] * b[2*i-3] * b[2*i+1], b[2*i-1] * ( a[2*i+1] + b[2*i] * b[2*i+1] ) + a[2*i] * b[2*i+1] ] ];
          i := i + 1;
        end do;
      end if;

      # then we create the part with general elements

      t := tformula['general'];  nopst := nops(t);
      a := [ seq( t[j][1], j=1..nopst ) ];
      b := [ seq( t[j][2], j=1..nopst ) ];

      j := 0;
      i := newformula['index'];
      newformula['general'] := [];

      if assigned( tformula['begin'] ) and is ( nops(tformula['begin']), odd ) then
        j := j + 1;
      end if;

      for k from 1 by 1 to `if`( type( nopst, even ), nopst/2, nopst ) do
        newformula['general'] := [ op(newformula['general']), [ - subs(i=2*i-1, a[(j+2) mod nopst + 1]) * subs(i=2*i, a[(j+3) mod nopst + 1]) * subs(i=2*i-3, b[j mod nopst + 1]) * subs(i=2*i+1, b[(j+4) mod nopst + 1]), subs(i=2*i-1, b[(j+2) mod nopst + 1]) * ( subs(i=2*i+1, a[(j+4) mod nopst + 1]) + subs(i=2*i, b[(j+3) mod nopst + 1]) * subs(i=2*i+1, b[(j+4) mod nopst + 1]) ) + subs(i=2*i, a[(j+3) mod nopst + 1]) * subs(i=2*i+1, b[(j+4) mod nopst + 1]) ] ];
        j := (j+1) mod nopst + 1;
      end do;

      if assigned( tformula['label'] ) then
        newformula['comment'] := cat( "odd contraction of continued fraction formula with label \"", tformula['label'], "\"" );
    end if;

    ##### ----- simregular -----

    elif form = 'simregular' then

      # first we check whether the formula is already in simregular form or not;
      # if this is the case, we just return the formula

      if type( tformula, 'simregular' ) then
        return ':-CFSF:-formula'( copy(tformula) );
      end if;

      # otherwise, we construct the front and begin parts

      t := [ seq( nthpelement( tformula, i ), i=0..nops(tformula['begin'])+1 ) ];
      a := [ seq( t[j][1], j=2..nops(t) ) ];
      b := [ seq( t[j][2], j=2..nops(t) ) ];

      if assigned( tformula['front'] ) then
        newformula['front'] := tformula['front'];
      end if;

      if assigned( tformula['factor'] ) then
        newformula['factor'] := tformula['factor'];
        newformula['begin'] := [ [ a[1] / ( b[1] * newformula['factor'] ), 1 ] ];
      else
        newformula['begin'] := [ [ a[1] / b[1], 1 ] ];
      end if;

      if assigned( tformula['begin'] ) then
        i := 2;
        while ( i-1 <= nops(tformula['begin']) ) do
          newformula['begin'] := [ op(newformula['begin']), [ a[i] / ( b[i-1] * b[i] ), 1 ] ];
          i := i + 1;
        end do;
      end if;

      # then we create the part with general elements

      t := tformula['general'];  nopst := nops(t);
      a := [ seq( t[j][1], j=1..nopst ) ];
      b := [ seq( t[j][2], j=1..nopst ) ];

      j := 0;
      i := newformula['index'];
      newformula['general'] := [];

      for k from 1 by 1 to nopst do
        newformula['general'] := [ op(newformula['general']), [ a[(j+1) mod nopst + 1] / ( subs(i=i-1, b[j mod nopst + 1]) * b[(j+1) mod nopst + 1] ) , 1 ] ];
        j := j mod nopst + 1;
      end do;

      if assigned( tformula['label'] ) then
        newformula['comment'] := cat( "simregular form of continued fraction formula with label \"", tformula['label'], "\"" );
      end if;

    ##### ----- Euler -----

    elif form = 'Euler' then

      newformula['type'] := 'series';
      newformula['index'] := ':-k';

      newformula['front'] := nthpelement( tformula, 0 );
      newformula['general'] := [ (-1)^(newformula['index']-1)
        * product( ':-CFSF:-nthpnumer'( aformula, ':-j' ), ':-j'=1..newformula['index'] )
        / ( ':-CFSF:-nthdenom'( aformula, newformula['index'] ) * ':-CFSF:-nthdenom'( aformula, newformula['index']-1 ) ) ];

      if assigned( tformula['label'] ) then
        newformula['comment' ] := cat( "Euler series of continued fraction formula with label \"", tformula['label'], "\"" );
      end if;

    else
      ERROR( "`%1' is not a valid transformation for formulas of type `%2'", form, tformula['type'] );
    end if;

    # before we return the newly created formula, we check whether some
    # of the begin elements which are conform with the general elements;
    # if so, we drop them from the formula

    if assigned( newformula['begin'] ) and nops( newformula['general'] ) = 1 then
      for i from nops( newformula['begin'] ) by -1 to 1 do
        if verify( subs( newformula['index']=i, newformula['general'][1] ), newformula['begin'][i] , 'list'('simplify') ) then
          newformula['begin'] := [ op(1..i-1,newformula['begin']) ];
        else
          break;
        end if;
      end do;
    end if;

    # we also try to simplify the front, begin and general elements

    try
      if assigned( newformula['front'] ) then newformula['front'] := simplify(newformula['front']); end if;
      if assigned( newformula['begin'] ) then newformula['begin'] := simplify(newformula['begin']); end if;
      if assigned( newformula['general'] ) then newformula['general'] := simplify(newformula['general']); end if;
    catch:
    end try;

    # finally, we do a safety check to make sure the formula validates
    # correctly; if not, we have an implementation problem (call it a bug)

    if not validate( newformula ) then
      ERROR( "validation failure, unable to create the transformation successfully" );
    end if;

    return ':-CFSF:-formula'( copy(eval(newformula)) );

  end proc;

end module;

##### ===== END OF CONTFRAC SUBMODULE =====
