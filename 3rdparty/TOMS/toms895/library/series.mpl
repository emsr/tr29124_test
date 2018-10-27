##### ============================
##### ===== SERIES SUBMODULE =====
##### ============================

series := module()

  description "CFSF series submodule";

  export
    argtypes, argnames, reqargs,
    create, formula,
    nthpelement,
    nthapprox,
    nthtail,
    transform;

  local
    setup, cleanup;

  option
    load = setup, unload = cleanup;

##### ----- LOADING AND UNLOADING -----

  setup := proc()
    description "setup procedure for CFSF[series]";
    USERINFO( 5, ':-CFSF', `setup procedure for CFSF[series]` );
    NULL
  end proc;

  cleanup := proc()
    description "cleanup procedure for CFSF[series]";
    USERINFO( 5, ':-CFSF', `cleanup procedure for CFSF[series]` );
    NULL
  end proc;

##### ----- CREATING AND VALIDATING FORMULAS -----

  argtypes := table([ op(op(op(cargtypes))),
    'even' = exprtype,
    'odd' = exprtype,
    'general' = list(exprtype)
  ]);

  argnames := map( op, { indices(argtypes) } );

  reqargs := { op(reqcargs), 'general' };

  create := proc( newformula::table )

    if not assigned( newformula['index'] ) then
      newformula['index'] := 'k';
    end if;

  end proc;

##### ----- RETRIEVING FORMULAS -----

  formula := proc( tformula::table )

    if nargs >= 2 and not type( args[2], set ) then
      ERROR( "`formula' expects its 2nd argument to be of type `set', but got `%1'", args[2] );
    end if;

    if nargs > 2 then
      ERROR( "`formula' expects 1 or 2 arguments, but got %1", nargs );
    end if;

  end proc;

##### ----- RETRIEVING PARTS OF FORMULAS -----

  nthpelement := proc( tformula::table, n::nonnegint )

    local tfactor;

    tfactor := `if`( assigned( tformula['factor'] ), tformula['factor'], 1 );

    if n = 0 then
      if assigned( tformula['front'] ) then
        eval( tformula['front'] );
      else
        eval( subs( tformula['index'] = 0, tformula['general'][1] ) * tfactor );
      end if;
    elif assigned( tformula['front'] ) then
      eval( subs( tformula['index'] = n, tformula['general'][(n-1) mod nops(tformula['general']) + 1] ) * tfactor );
    else
      eval( subs( tformula['index'] = n, tformula['general'][n mod nops(tformula['general']) + 1] ) * tfactor );
    end if;

  end proc;

  ##### ----- EVALUATING FORMULAS -----

  nthapprox := proc( tformula::table, n::nonnegint, tsubs::table )

    local generallength, tsubstitutions, expr, t, i;

    generallength := nops( tformula['general'] );
    tsubstitutions := constructsubs( tformula, tsubs );

    # using backward recurrence (not using nthpelement, since
    # we do not want to multiply every term with the factor)

    # NOTE : we do an extra eval since the action of substitution
    # is not followed by evaluation (according to ?subs)

    expr := 0;

    t := `if`( assigned( tformula['front'] ), 1, 0 );

    for i from n by -1 to t do
      expr := eval( subs( tsubstitutions union { tformula['index'] = i }, tformula['general'][(i-t) mod generallength + 1] ) ) + expr;
      if assigned( tsubs['variable'] ) and type( subs['variable'], complex(float) ) then
        expr := evalf( expr );
      elif defaults[evaluate_simplify] > 0 and i mod defaults[evaluate_simplify] = 0 then
        expr := simplify( expr );
      end if;
    end do;

    if assigned( tformula['factor'] ) then
      expr := eval( subs( tsubstitutions, tformula['factor'] ) ) * expr;
    end if;

    if assigned( tformula['front'] ) then
      expr := eval( subs( tsubstitutions, tformula['front'] ) ) + expr;
    end if;

    # if the variable is of type float, then we do an extra evalf
    # to eliminate constants from the output (e.g. Pi)

    if assigned( tsubs['variable'] ) and type( tsubs['variable'], complex(float) ) then
      expr := evalf( expr );
    end if;

    # finally, we normalize and sort the output

    return sort( simplify( expr ) );

  end proc;

##### ----- WORKING WITH TAILS -----

  nthtail := proc( tformula::table, n::nonnegint )

    local newformula, t;

    # first we copy some arguments from the original formula

    newformula := table();
    newformula['type'] := tformula['type'];
    newformula['index'] := tformula['index'];
    newformula['variable'] := tformula['variable'];

    if assigned( tformula['parameters'] ) then
      newformula['parameters'] := tformula['parameters'];
    end if;

    # 'front' never reappears and 'factor' always stays

    if assigned( tformula['factor'] ) then
      newformula['factor'] := tformula['factor'];
    end if;

    # now we adjust the 'general' argument
    # depending on the value of n

    t := `if`( assigned( tformula['front'] ), 0, 1 );
    t := (n+t) mod nops( tformula['general'] ) + 1;

    if t = 1 then
      newformula['general'] := tformula['general'];
    else
      newformula['general'] := [ op(tformula['general'])[t..-1], op(tformula['general'])[1..t-1] ];
    end if;

    newformula['general'] := subs( newformula['index'] = newformula['index'] + n+1, newformula['general'] );

    # we do a safety check to make sure the formula validates correctly;
    # if not, we have found a bug in the implementation

    if not validate( newformula ) then
      ERROR( "validation failure, unable to create the tail successfully" );
    end if;

    if assigned( tformula['label'] ) then
      newformula['comment'] := cat( convert( n, string ), "-th tail of series formula with label \"", tformula['label'], "\"" );
    end if;

    # finally, wrap the table in a 'formula' expression

    return ':-CFSF:-formula'( copy(newformula) );

  end proc;

##### ----- TRANSFORMATIONS -----

  transform := proc( aformula::uneval, form::name )

    local tformula, newformula, t, nopst, i, j, k;

    tformula := eval( aformula );

    # first we copy some arguments from the original formula

    newformula := table();
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

    ##### ----- Euler -----

    if form = 'Euler' then

      newformula['type'] := 'contfrac';

      # first we construct the front and begin parts

      newformula['front'] := nthpelement( tformula, 0 );
      newformula['begin'] := [ [ nthpelement( tformula, 1 ), 1 ] ];

      # then we create the part with general elements

      t := tformula['general']; nopst := nops(t);

      j := 0;
      i := newformula['index'];
      newformula['general'] := [];

      if not assigned( tformula['front'] ) then
        j := j+1;
      end if;

      for k from 1 by 1 to nopst do
        newformula['general'] := [ op(newformula['general']), [ - t[(j+1) mod nopst + 1] / subs(i=i-1, t[j mod nopst + 1]), 1 + t[(j+1) mod nopst + 1] / subs(i=i-1, t[j mod nopst + 1]) ] ];
        j := j mod nopst + 1;
      end do;

      if assigned( tformula['label'] ) then
        newformula['comment'] := cat( "Euler form of series formula with label \"", tformula['label'], "\"" );
    end if;

    else
      ERROR( "`%1' is not a valid transformation for formulas of type `%2'", form, tformula['type'] );
    end if;

    # finally, we do a safety check to make sure the formula validates
    # correctly; if not, then there is a bug in this routine

    if not validate( newformula ) then
      ERROR( "validation failure, unable to create the transformation successfully" );
    end if;

    return ':-CFSF:-formula'( copy(simplify(eval(newformula))) );

  end proc;

end module;

##### ===== END OF SERIES SUBMODULE =====
