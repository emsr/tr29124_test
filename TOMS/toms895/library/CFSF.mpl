# to (re)create the CFSF package you have to execute 'maple -s -q CFSF.mpl'
# on the command line and then you will get two files CFSF.lib and CFSF.ind,
# next you have to add the path of these files to the variable 'libname',
# or you may add the line 'libname := libname, "path-to-CFSF.lib"' to the
# initialization files (.mapleinit for unix or maple.ini for windows)

# NOTE : for more technical information about the structure and/or layout
# of this package, take a look at the Shapes example package provided in
# Maple by the shapes.mpl file ...

##### ----- LIBRARY CREATION SETTINGS -----

# the following disables the printing of warning messages,
# it is used to avoid warnings about implicitly declared
# local formulas when creating the library - NOTE : these
# local formulas no longer exist, so we enable it again

# interface( warnlevel = 0 ):

# the infolevel mechanism is used to show some information
# about the declared formulas when creating the library

infolevel[':-CFSF'] := 3:

# the following alias for userinfo is used only for typography purposes,
# i.e. to make it look more like WARNING and ERROR (which for the same
# reason is used instead of error)

alias( USERINFO = userinfo ):

##### ============================
##### ===== CFSF MAIN MODULE =====
##### ============================

module CFSF()

  description "CFSF package";

  export
    functions,
    create, add_formula, remove_formula,
    query, formula,
    nthpelement, nthpnumer, nthpdenom,
    nthnumer, nthdenom, nthapprox,
    nthtail, tailestimate,
    transform, varsubs, equivalent,
    default, build, init;

  local
    CFSF, setup, cleanup,
    zerotype, exprtype, symboltype, functiontype,
    cargtypes, reqcargs, validate,
    processargs, constructsubs, dispatch,
    bname, bdate, bplatform, bversion, btarget,
    formulas, defaults;

  option
    package, load = setup, unload = cleanup;

##### ----- LOADING AND UNLOADING -----

  setup := proc()

    description "setup procedure for CFSF";

    global `type/formula`, `print/formula`;

    USERINFO( 5, ':-CFSF', `setup procedure for CFSF` );

    # here we define the formula type using the classical extension mechanism,
    # it uses the 'validate' procedure for type checking

    `type/formula` := proc( aformula )

      local tformula;

      # USERINFO( 5, ':-CFSF', `called with arguments`, args );

      tformula := eval(aformula);

      # not type( aformula, table )
      # and ( type( aformula, symbol ) or type( aformula, function ) )
      # and
      type( tformula, table )
      and validate( tformula, 'warnings', 'halt' );

    end proc;

    # the following extension to the print routine allows us to control
    # the output when working with expressions of the form 'formula'(...),
    # so we can hide the table structure that is used underneath
    # (which is useful for readability and the like)

    `print/formula` := proc( aformula )

      local tformula;

      # USERINFO( 5, ':-CFSF', `called with arguments`, args );

      # NOTE : we only want pretty printing when the CFSF package is loaded
      # this is needed because of the gc() problem (see cleanup)

      if member( ':-CFSF', packages() ) then

        if ( ( type( aformula, function ) and evalb( op(0,aformula) = TABLE ) ) or type( aformula, table ) )
        and ( type( op(aformula), list ) or type( op(op(aformula)), list ) ) then

          # note that we have to explicitly recreate the table structure since
          # aformula is of the form TABLE([...]) and TABLE is an undocumented
          # name that is reserved for internal use (see ?TABLE) - in case of
          # a regular table structure, we can use its structure immediately
          # (for internal purposes only)

          tformula := `if`( type( aformula, table ), aformula, table( op(aformula) ) );

          # NOTE : checking for a formula from the formulas table is done
          # by comparing all arguments (i.e. the complete table), which
          # may be an expensive operation

          if validate( tformula )
          then
            # if infolevel[':-CFSF'] = 1 or infolevel[':-CFSF'] = 2 then
            #   ':-CFSF:-formula'( op(aformula) );
            if assigned( tformula['label'] )
            and assigned( formulas[tformula['label']] )
            and verify( tformula, formulas[tformula['label']], 'table' )
            then
              ':-CFSF:-formula'( tformula['label'] );
            else
              ':-CFSF:-formula'([ 'type'=tformula['type'], `...` ]);
              # ':-CFSF:-formula'( op(op(tformula)) );
            end if;
          else
            WARNING( "this is not a valid CFSF formula" );
            ':-formula'( op(op(tformula)) );
          end if;

        elif type( aformula, string ) and assigned( formulas[aformula] ) then

          ':-CFSF:-formula'( aformula );

        else

	      # NOTE : if the 'formula' procedure name was used to create a named
          # object before the CFSF package was loaded, it cannot be reused
          # afterwards and has to be explicitly recreated (that is, reassigned)

          USERINFO( 4, ':-CFSF', `this is not a valid CFSF formula` );
          WARNING( "this is not a valid CFSF formula" );
        end if;

      # else
      #   ':-formula'( aformula );
      end if;

    end proc;

    NULL

  end proc;

  # NOTE : when garbage collection is done (either automatically as specified
  # by kernelopts(gcfreq), or by manually executing gc()) the cleanup routine
  # is called, even when we are still using the package - until we figure out
  # why this happens (and fix it), nothing should actually be cleaned up here
  # to avoid possible problems

  cleanup := proc()

    description "cleanup procedure for CFSF";

    global `type/formula`, `print/formula`;

    USERINFO( 5, ':-CFSF', `cleanup procedure for CFSF` );

    # `type/formula` := evaln( `type/formula` );
    # `print/formula` := evaln( `print/formula` );

    NULL

  end proc;

##### ----- EXPRTYPE AND ARGNAMES -----

  zerotype := Or( poszero, negzero, cx_zero );
  exprtype := { `+`, `*`, `^`, function, symbol, complex(numeric) };
  symboltype := And( atomic, symbol, Not( constant ) );
  functiontype := { constant, procedure };

  # cargtypes defines the types for the common arguments of the formulas;
  # it is extended in the submodules and used by the 'validate' procedures

  cargtypes := table([
    'label' = string, 'booklabel' = string,
    'front' = exprtype, 'factor' = exprtype,
    'variable' = symboltype, 'index' = symboltype,
    'parameters' = set(symboltype),
    'function' = { functiontype, set( functiontype ) }, 'lhs' = exprtype,
    'category' = string, 'comment' = string,
    'constraints' = set( { exprtype::property, relation, logical } )
  ]);

  # reqcargs defines the arguments that are absolutely required and therefore
  # must be defined (either explicitly, otherwise implicitly by providing a
  # default value for them in the 'create' routines of the submodules)
  # (it can also include argnames which are not yet defined here,
  # in this case they must be explicitly defined in every submodule)

  reqcargs := { 'variable', 'index' };

##### ----- VALIDATING FORMULAS -----

  # the 'validate' procedure provides strong validation rules for formulas:
  # it checks for missing arguments, wrong argument types, and/or assignments
  # of arguments; it is used by 'type/formula' (so wherever ::formula occurs),
  # 'print/formula' (so every time a formula is printed) as well as
  # by 'create' and 'formula' (and possibly other procedures)
  # NOTE : using this procedure implies safety at the cost of frequently
  # and repeatedly occuring type checking

  validate := proc( tformula::table )

    local check, alabel, labels, i;

    USERINFO( 5, ':-CFSF', `validation procedure for CFSF formulas` );

    # first we check whether it is a valid CFSF formula

    if not assigned( tformula['type'] )
    or not member( tformula['type'], CFSF )
    then
      USERINFO( 4, ':-CFSF', `this is not a valid CFSF formula` );
      return false;
    end if;

    # then we check for missing required arguments; if 'missing' is specified
    # as the second argument, then we return the list of missing arguments
    # (if any, otherwise we return the empty list)

    labels := map( op, { indices(tformula) } );
    check := tformula['type']:-reqargs minus labels;

    if nargs = 2 and args[2] = 'missing' then
      return check;
    elif nops( check ) > 0 then
      USERINFO( 4, ':-CFSF', `missing required arguments`, check, `for formulas of type`, tformula['type'] );
      USERINFO( 4, ':-CFSF', `this is not a valid CFSF formula` );
      return false;
    end if;

    # if 'warnings' is specified as the second argument, we issue a warning
    # when either the variable, index or any of the parameters has associated
    # assumptions or a value assigned to it, or when the factor is zero;
    # if furthermore 'halt' is specified as a third argument, we stop here

    check := true;

    if not evalb( tformula['variable'] in {unames()} ) then
      if hasassumptions( eval(tformula['variable']) ) then
        USERINFO( 4, ':-CFSF', `the variable used in this formula has associated assumptions` );
        if nargs >= 2 and args[2] = 'warnings' then
          WARNING( "the variable used in this formula has associated assumptions" );
        end if;
      else # elif assigned( eval(tformula['variable']) ) then
        USERINFO( 4, ':-CFSF', `the variable used in this formula has a value assigned to it` );
        if nargs >= 2 and args[2] = 'warnings' then
          WARNING( "the variable used in this formula has a value assigned to it" );
        end if;
      end if;
      check := false;
    end if;

    if not evalb( tformula['index'] in {unames()} ) then
      if hasassumptions( eval(tformula['index']) ) then
        USERINFO( 4, ':-CFSF', `the index used in this formula has associated assumptions` );
        if nargs >= 2 and args[2] = 'warnings' then
          WARNING( "the index used in this formula has associated assumptions" );
        end if;
      else # elif assigned( eval(tformula['index']) ) then
        USERINFO( 4, ':-CFSF', `the index used in this formula has a value assigned to it` );
        if nargs >= 2 and args[2] = 'warnings' then
          WARNING( "the index used in this formula has a value assigned to it" );
        end if;
      end if;
      check := false;
    end if;

    if assigned( tformula['parameters'] ) then
      for i from 1 to nops( tformula['parameters'] ) do
        if not evalb( tformula['parameters'][i] in {unames()} ) then
          if hasassumptions( eval(tformula['parameters'][i]) ) then
            USERINFO( 4, ':-CFSF', `one of the parameters used in this formula has associated assumptions` );
            if nargs >= 2 and args[2] = 'warnings' then
              WARNING( "one of the parameters used in this formula has associated assumptions" );
            end if;
          else
            USERINFO( 4, ':-CFSF', `one of the parameters used in this formula a value assigned to it` );
            if nargs >= 2 and args[2] = 'warnings' then
              WARNING( "one of the parameters used in this formula has a value assigned to it" );
            end if;
          end if;
          check := false;
        end if;
      end do;
    end if;

    if assigned( tformula['factor'] ) then
      if type( tformula['factor'], zerotype ) then
        USERINFO( 4, ':-CFSF', `the factor used in this formula is zero` );
        if nargs >= 2 and args[2] = 'warnings' then
          WARNING( "the factor used in this formula is zero" );
        end if;
          check := false;
      end if;
    end if;

    if nargs = 3 and args[2] = 'warnings' and args[3] = 'halt' then
      if not check then
        USERINFO( 4, ':-CFSF', `this is not a valid CFSF formula` );
      end if;
      return check;
    end if;

    # now check whether all arguments are of the correct type
    # NOTE : we use eval(,1) to avoid returning false in case
    # of associated assumptions and/or assigned values

    for alabel in ( labels minus { 'type' } ) do
      if alabel in tformula['type']:-argnames then
        check := type( eval(tformula[alabel],1), tformula['type']:-argtypes[alabel] );
        if check = false then
          USERINFO( 4, ':-CFSF', `expecting an argument of type`, tformula['type']:-argtypes[alabel], `for argument`, alabel, `of formula type`, tformula['type'], `but got`, tformula[alabel] );
          USERINFO( 4, ':-CFSF', `this is not a valid CFSF formula` );
          return false;
        end if;
      else
        USERINFO( 4, ':-CFSF', alabel, `is not a valid argument for formulas of type`, tformula['type'] );
        USERINFO( 4, ':-CFSF', `this is not a valid CFSF formula` );
        return false;
      end if;
    end do;

    # for the constraints, we need extra checking
    # NOTE : currently, we only issue a warning

    if assigned( tformula['constraints'] ) then
      for i in tformula['constraints'] do
        if ( type( i, function ) and op(0,i) = `in` and not type( op(2,i), property ) )
        or ( type( i, `::` ) and not type( op(2,i), property ) ) then
          USERINFO( 4, ':-CFSF', i, `is not a valid constraint` );
          if nargs >= 2 and args[2] = 'warnings' then
            WARNING( "one of the constraints from this formula is not valid" );
          end if;
          break;
        end if;
      end do;
    end if;

    # we call the formula-specific 'validate' procedure, if it exists

    if member( ':-validate', tformula['type'] ) then
      check := tformula['type']:-validate( tformula );
    end if;

    if check then
      # USERINFO( 4, ':-CFSF', `this is a valid CFSF formula` );
    else
      USERINFO( 4, ':-CFSF', `this is not a valid CFSF formula` );
    end if;

    return check;

  end proc;

##### ----- SUBMODULES -----

  # the following submodule is used as a wrapper around the submodules
  # for each of the supported formula types, which makes it possible
  # to use Maple-like names for these submodules and avoids having
  # to export them (i.e. they will not be reinterpreted)

  CFSF := module()

    export
      series, contfrac;

$include "series.mpl"
$include "contfrac.mpl"

  end module;

  # also include the (new) functions and constants submodule

$include "functions.mpl"

##### ----- CREATING FORMULAS -----

  create := proc( newftype::symbol )

    description "procedure for creating formulas";

    local ftype, newformula, i;

    # first, check whether it is a known formula type

    if not member( newftype, CFSF, 'ftype' ) then
      ERROR( "`%1' is not a valid CFSF formula type", aformula );
    else
      USERINFO( 4, ':-CFSF', `the formula is of type`, ftype );
    end if;

    # next, process the supplied arguments

    newformula := table();
    newformula['type'] := ftype;

    for i in args[2..nargs] do
      if not type( i, equation ) then
        ERROR( "`create' expects its arguments to be of the form `<argument> = <expression>', but got `%1'", i );
      elif not member( op(1,i), CFSF[ftype]:-argnames ) then
        ERROR( "`%1' is not a valid argument for formulas of type `%2'", op(1,i), ftype );
      elif not type( op(2,i), CFSF[ftype]:-argtypes[op(1,i)] ) then
        ERROR( "`create' expects an expression of type `%1' for argument `%2', but got `%3'", CFSF[ftype]:-argtypes[op(1,i)], op(1,i), op(2,i) );
      else
        newformula[op(1,i)] := op(2,i);
        USERINFO( 4, ':-CFSF', `the formula has argument`, op(1,i), `with value`, op(2,i) );
      end if;
    end do;

    # implicitly define values for the required arguments (from reqcargs)

    if not assigned( newformula['variable'] ) then
      newformula['variable'] := 'z';
    end if;

    # next, we call the 'create' procedure of the specified formula type,
    # which allows us to do some formula-specific things (such as setting
    # the required arguments if they are not specified explicitly),
    # if it exists

    if member( ':-create', CFSF[ftype] ) then
      CFSF[ftype]:-create( newformula );
    end if;

    # now, check for construction errors and/or transform the arguments
    # when necessary/possible (e.g., create a `general' argument out of
    # `even' and `odd' arguments, ...)

    if not ( assigned( newformula['even'] ) or assigned( newformula['odd'] ) or assigned( newformula['general'] ) )
    or not ( ( assigned( newformula['even'] ) and assigned( newformula['odd'] ) ) or assigned( newformula['general'] ) )
    or ( ( assigned( newformula['even'] ) or assigned( newformula['odd'] ) ) and assigned( newformula['general'] ) )
    then
      ERROR( "a formula of type `%1' must have either both 'even' and 'odd' arguments or a 'general' argument", newformula['type'] );
    elif ( assigned( newformula['even'] ) and assigned( newformula['odd'] ) ) then
      USERINFO( 4, ':-CFSF', `replacing 'even' and 'odd' by 'general'` );
      if not assigned( newformula['begin'] ) or type( nops( newformula['begin'] ), even ) then
        newformula['general'] := [ newformula['odd'], newformula['even'] ];
      else
        newformula['general'] := [ newformula['even'], newformula['odd'] ];
      end if;
      unassign( 'newformula[even]' );
      unassign( 'newformula[odd]' );
    end if;

    # the 'validate' procedure is called to check whether we are not
    # missing any required arguments (no need to recheck everything)
    # NOTE : this should be an unnecessary safety check, since
    # everything should be covered already in the above part

    i := validate( newformula, 'missing' );
    if nops(i) > 0 then
      ERROR( "missing required arguments `%1' for formulas of type `%2'", i, newformula['type'] );
    end if;

    # EXTRA : for safety reasons, we do not allow a formula to be created when
    # either the variable, the index or the parameters (if any) has associated
    # assumptions or an assigned value, or when the factor is zero

    if validate( newformula, 'warnings', 'halt' ) = false then
      # ERROR( "the variable, the index and all parameters must not have any associated assumptions or assigned values, and the factor must be nonzero" );
      ERROR( "warnings encountered, unable to create the formula successfully" );
    end if;

    # finally, if a label is given, we add it to the formulas table
    # if the variable 'defaults[add_formulas]' has been set to 'true'

    if assigned( newformula['label'] ) then
      if defaults['add_formulas'] = false then
        USERINFO( 3, ':-CFSF', `the formula has not been added because it is not allowed` );
      elif assigned( formulas[newformula['label']] ) then
        USERINFO( 3, ':-CFSF', `the formula could not be added because there already exists a formula with label`, newformula['label'] );
      else
        formulas[ newformula['label'] ] := newformula;
        USERINFO( 3, ':-CFSF', `the formula has been added and is accessible as`, convert( cat( "formula(\"", newformula['label'], "\")" ), symbol ) );
        # return ':-CFSF:-formula'( newformula['label'] );
      end if;
    end if;

    # we return the table wrapped in a 'formula' expression so that
    # the `print/formula' procedure will be called to format the output

    return ':-CFSF:-formula'( copy(newformula) );

  end proc;

##### ----- TABLE OF KNOWN FORMULAS -----

  # when you create a formula that contains a label, it will
  # be added to the formulas table, provided that both the variable
  # 'default(add_formulas)' has been set to 'true' and also that
  # this label does not already appear in the formulas table

  formulas := table();

  # the initial settings for the formulas table is located in a separate file

$include "formulas.mpl"

  # the following functionality can be used to add and/or remove
  # formulas from the formulas table

  add_formula := proc( aformula::formula, alabel::string )

    local tformula;

    tformula := copy(aformula);

    if defaults['add_formulas'] = false then
      USERINFO( 3, ':-CFSF', `the formula has not been added because it is not allowed` );
      WARNING( "the formula has not been added because it is not allowed" );
    elif assigned( formulas[alabel] ) then
      USERINFO( 3, ':-CFSF', `the formula could not be added because there already exists a formula with label \"%1\"`, alabel );
      WARNING( "the formula could not be added because there already exists a formula with label \"%1\"", alabel );
    else
      if assigned( tformula['label'] ) and not tformula['label'] = alabel then
        # WARNING( "the formula previously had label \"%1\"", tformula['label'] );
      end if;
      tformula['label'] := alabel;
      formulas[alabel] := tformula;
      USERINFO( 3, ':-CFSF', `the formula has been added and is accessible as`, convert( cat( "formula(\"", alabel, "\")" ), symbol ) );
    end if;

  end proc;

  remove_formula := proc( alabel::string )

    if assigned( formulas[alabel] ) then
      if defaults['add_formulas'] = true then
        formulas[alabel] := 'formulas[alabel]';
        USERINFO( 3, `the formula with label`, alabel, `has been removed` );
      else
        WARNING( "the formula has not been removed since it is not allowed" );
      end if;
    else
      WARNING( "no formula with label \"%1\" exists", alabel );
    end if;

  end proc;



##### ----- QUERYING THE FORMULAS TABLE -----

  query := proc()

    local function, category, flabels, selected, i;

    flabels := map( op, { indices(formulas) } );

    USERINFO( 4, ':-CFSF', `starting with`, nops( labels ), `formulas` );

    if nargs = 0 then
      return flabels;
    elif nargs > 2 then
      ERROR( "`query' expects up to 2 arguments, but got %1", nargs );
    end if;

    if type( args[1], equation ) then
      if op(1,args[1]) = ':-function' and type( op(2,args[1]), cargtypes[':-function'] ) then
        function := op(2,args[1]);
      elif op(1,args[1]) = ':-category' and type( op(2,args[1]), cargtypes[':-category'] ) then
        category := op(2,args[1]);
      else
        ERROR( "`query' expects its first argument to be of the form `'function' = <symbol>' or `'category' = <string>', but got `%1'", args[1] );
      end if;
    else
      ERROR( "`query' expects its first argument to be of the form `'function' = <symbol>' or `'category' = <string>', but got `%1'", args[1] );
    end if;

    if nargs = 2 then
      if type( args[2], equation ) then
        if op(1,args[2]) = ':-function' and type( op(2,args[2]), cargtypes[':-function'] ) then
          if not assigned( function ) then
            function := op(2,args[2]);
          else
            ERROR( "`query' allows only one argument of the form `'function' = <symbol>" );
          end if;
        elif op(1,args[2]) = ':-category' and type( op(2,args[2]), cargtypes[':-category'] ) then
          if not assigned( category ) then
            category := op(2,args[2]);
          else
            ERROR( "`query' allows only one argument of the form `'category' = <string>" );
          end if;
        else
        ERROR( "`query expects its second argument to be of the form `'function' = <symbol>' or `'category' = <string>', but got `%1'", args[2] );
        end if;
      else
        ERROR( "`query expects its second argument to be of the form `'function' = <symbol>' or `'category' = <string>', but got `%1'", args[2] );
      end if;
    end if;

    if assigned( function ) then
      selected := {};
      for i from 1 to nops( flabels ) do
        if assigned( formulas[flabels[i]][':-function'] )
        and ( ( type( formulas[flabels[i]][':-function'], set ) and function in formulas[flabels[i]][':-function'] )
        or ( formulas[flabels[i]][':-function'] = function ) )
	then
          selected := selected union {flabels[i]};
        end if;
      end do;
      flabels := selected;
      USERINFO( 4, ':-CFSF', `found`, nops( selected ), `formula(s) with function =`, function );
    else
      USERINFO( 4, ':-CFSF', `no function query requested` );
    end if;

    if assigned( category ) then
      selected := {};
      for i from 1 to nops( flabels ) do
        if assigned( formulas[flabels[i]][':-category'] )
        and searchtext( category, formulas[flabels[i]][':-category'] ) <> 0
        then
          selected := selected union {flabels[i]};
        end if;
      end do;
      flabels := selected;
      USERINFO( 4, ':-CFSF', `found`, nops( selected ), `formula(s) with category =`, category );
    else
      USERINFO( 4, ':-CFSF', `no category query requested` );
    end if;

    return flabels;

  end proc;

##### ----- RETRIEVING FORMULAS -----

  formula := proc( alabel )

    description "function to get a formula from the formulas table";

    local tformula, tparams, tconstraints, amodification, i, j;

    # USERINFO( 5, ':-CFSF', `in formula procedure, called with arguments`, args );

    if nargs < 1 or nargs > 3 then
      ERROR( "`formula' expects 1 to 3 arguments, but got %1", nargs );
    end if;

    # if the first argument is of type table, the formula was already
    # constructed, so we just return the table it contains (this way,
    # one can see the underlying table structure using eval)

    if type( alabel, 'table' ) then
      # USERINFO( 5, ':-CFSF', `in formula procedure, returning table` );
      return eval(alabel);
    end if;

    # otherwise the user wants to retrieve a formula from the formulas
    # table; in this case we first create a copy of the requested formula,
    # since we do not want people to work with the original one

    if not type( alabel, cargtypes['label'] ) then
      ERROR( "`formula' expects its 1st argument to be of type `%1', but got `%2'", cargtypes['label'], alabel );
    elif not assigned( formulas[alabel] ) then
      ERROR( "no formula with label \"%1\" exists", alabel );
    else
      tformula := copy( eval(formulas[alabel]) );
    end if;

    # we call the formula-specific 'formula' procedure, if it exists

    if member( ':-formula', tformula['type'] ) then
      tformula['type']:-formula( tformula, args[2..nargs] );
    end if;

    # we issue a warning if the formula variable or its parameters
    # have a value assigned to it (use at your own risk)

    validate( tformula, 'warnings' );

    # if a set of equations is given for (some of) the parameters,
    # we substitute each given parameter (the lhs of the equation)
    # with its supplied value (the rhs of the equation)

    if nargs >= 2 and type( args[2], set ) then

      if assigned( tformula['parameters'] ) and nops( args[2] ) > 0 then

        tparams := args[2];

        for i in tparams do
          if not type( i, equation ) then
            ERROR( "`formula' expects its 2nd argument to be of type `set' with members of the form `<parameter> = <expression>', but got `%1'", i );
          elif not member( op(1,i), tformula['parameters'] ) then
            tparams := tparams minus { i };
            WARNING( "`%1' is not a valid parameter for the requested formula", op(1,i) );
            WARNING( "ignoring substitution `%1'", i );
          else
            tformula['parameters'] := tformula['parameters'] minus { op(1,i) };
          end if;
        end do;

        # here, we do all the parameter substitutions simultaneously,
        # remove all substituted parameters, and check whether we need
        # to add any new parameters by looking at the right hand sides
        # of the substitutions (this allows renaming of parameters)

        if nops( tparams ) > 0 then

          if assigned( tformula['label'] ) then
            tformula['label'] := 'tformula[label]';
          end if;

          if assigned( tformula['booklabel'] ) then
            tformula['booklabel'] := 'tformula[booklabel]';
          end if;

          if assigned( tformula['constraints'] ) then
            tconstraints := tformula['constraints'];
          end if;

          tformula := subs( tparams, eval(tformula) );
          tformula['parameters'] := tformula['parameters'] union indets( map( rhs, tparams ), symboltype );
          tformula['parameters'] := tformula['parameters'] minus { tformula['variable'] };

          if nops( tformula['parameters'] ) = 0 then
            tformula['parameters'] := 'tformula[parameters]';
          end if;

          # when parameters have been substituted, check whether we
          # can drop a constraint that is always valid or whether we
          # should return an error because it evaluates to false

          if assigned( tformula['constraints'] ) then

            for i in tformula['constraints'] do
              if indets( i, symboltype ) = {} then
                if is( i ) = true then
                  tformula['constraints'] := tformula['constraints'] minus { i };
                  USERINFO( 4, ':-CFSF', `removing constraint`, i, `since it is always valid` );
                elif is( i ) = false then
                  ERROR( "constraint `%1' violated by parameter substitutions", i );
                end if;
              end if;
            end do;

            if nops( tformula['constraints'] ) = 0 then
              tformula['constraints'] := 'tformula[constraints]';
            end if;

          end if;

        end if;

      elif nops( args[2] ) > 0 then

        WARNING( "the requested formula has no parameters" );
        WARNING( "ignoring 2nd argument" );

      end if;

    end if;

    # again we wrap the table in a 'formula' expression to allow
    # custom output through `print/formula`

    return ':-CFSF:-formula'( copy(tformula) );

  end proc;

##### ----- ARGUMENT PROCESSING AND DISPATCHING -----

  # the procedure below is used by dispatch for processing supplementary
  # arguments when calling the nth* evaluation procedures; when all arguments
  # parse successfully, a table is returned with the substitutions that should
  # be executed during evaluation; otherwise an error is triggered

  processargs := proc( pname, fname, tformula )

    local tsubs, tparams, tsubstitutions, i;

    tsubs := table();

    if nargs = 3 then
      return eval(tsubs);
    end if;

    # check the given arguments for correctness

    if evalb( 'modification' in tformula['type']:-argnames ) then
      if nargs > 6 then
        ERROR( "`%1' expects 2 to 5 arguments, but got %2", fname, nargs-1 );
      end if;
    elif nargs > 5 then
      ERROR( "`%1' expects 2 to 4 arguments, but got %2", fname, nargs-1 );
    end if;

    # NOTE : the eval in comparison op(1,args[4]) = eval(tformula['variable'])
    # is needed because is(z~='z') evaluates to false when z has associated
    # assumptions assigned to it ...

    if nargs = 4 then

      if type( args[4], set ) then
        tparams := args[4];
      elif evalb( 'modification' in tformula['type']:-argnames ) and type( args[4], equation ) and op(1,args[4]) = 'modification' then
        tsubs['modification'] := op(2,args[4]);
      elif type( args[4], equation ) and op(1,args[4]) = eval(tformula['variable']) then
        tsubs['variable'] := op(2,args[4]);
      elif not evalb( 'modification' in tformula['type']:-argnames ) then
        ERROR( "`%1' expects its 3rd argument to be of type `set' or of the form `<variable> = <expression>', but got `%2'", fname, args[4] );
      else
        ERROR( "`%1` expects its 3rd argument to be of type `set` or of the form `'modification' = <expression>' or `<variable> = <expression>', but got `%2'", fname, args[4] );
      end if;

    elif nargs = 5 then

      if type( args[4], set ) then
        tparams := args[4];
      elif evalb( 'modification' in tformula['type']:-argnames ) and type( args[4], equation ) and op(1,args[4]) = 'modification' then
        tsubs['modification'] := op(2,args[4]);
      elif not evalb( 'modification' in tformula['type']:-argnames ) then
        ERROR( "`%1' expects its 3rd argument to be of type `set', but got `%2'", fname, args[4] );
      else
        ERROR( "`%1' expects its 3rd argument to be of type `set' or of the form `'modification' = <expression>', but got `%2'", fname, args[4] );
      end if;

      if evalb( 'modification' in tformula['type']:-argnames ) and type( args[5], equation ) and op(1,args[5]) = 'modification' and not assigned( tsubs['modification'] ) then
        tsubs['modification'] := op(2,args[5]);
      elif type( args[5], equation ) and op(1,args[5]) = eval(tformula['variable']) then
        tsubs['variable'] := op(2,args[5]);
      elif not evalb( 'modification' in tformula['type']:-argnames ) or ( evalb( 'modification' in tformula['type']:-argnames ) and assigned( tsubs['modification'] ) ) then
        ERROR( "`%1' expects its 4th argument to be of the form `<variable> = <expression>', but got `%2'", fname, args[5] );
      else
        ERROR( "`%1' expects its 4th argument to be of the form `'modification' = <expression>' or `<variable> = <expression>', but got `%1'", args[5] );
      end if;

    elif evalb( 'modification' in tformula['type']:-argnames ) and nargs = 6 then

      if type( args[4], set ) then
        tparams := args[4];
      else
        ERROR( "`%1' expects its 3rd argument to be of type `set', but got `%2'", fname, args[4] );
      end if;

      if type( args[5], equation ) and op(1,args[5]) = 'modification' then
        tsubs['modification'] := op(2,args[5]);
      else
        ERROR( "`%1' expects its 4th argument to be of the form `'modification' = <expression>', but got `%2'", fname, args[5] );
      end if;

      if type( args[6], equation ) and op(1,args[6]) = eval(tformula['variable']) then
        tsubs['variable'] := op(2,args[6]);
      else
        ERROR( "`%1' expects its 5th argument to be of the form `<variable> = <expression>', but got `%2'", fname, args[6] );
      end if;

    end if;

   # check the substitutions for the parameters

   if assigned( tparams ) and nops( tparams ) > 0 then

     if assigned( tformula['parameters'] ) then

       tsubs['parameters'] := {};

       for i from 1 to nops( tparams ) do
         if not type( op(i,tparams), equation ) then
           ERROR( "`%1' expects its 3rd argument to be of type `set' with members of the form `<parameter> = <expression>', but got `%2'", fname, op(i,tparams) );
         elif not member( op(1,op(i,tparams)), tformula['parameters'] ) then
           WARNING( "`%1' is not a valid parameter for this formula", op(1,op(i,tparams)) );
           WARNING( "ignoring substition `%1'", op(i,tparams) );
         else
           tsubs['parameters'] := tsubs['parameters'] union { op(i,tparams) };
         end if;
       end do;

     else
       WARNING( "the formula has no parameters" );
       WARNING( "ignoring 3rd argument" );
     end if;

   end if;

   # when a numeric value for the variable is given, check whether
   # the constraints are valid, and show a warning if it is not so

   if assigned( tsubs['variable'] ) and type( tsubs['variable'], complex(extended_numeric) ) and assigned( tformula['constraints'] ) then

     tsubstitutions := { tformula['variable'] = tsubs['variable'] };

     if assigned( tsubs['parameters'] ) then
       tsubstitutions := tsubstitutions union tsubs['parameters'];
     end if;

     for i in tformula['constraints'] do

     # NOTE : we use a try-catch block here to avoid possible
     # (and/or unknown) problems with checking the constraints
     # (the extra eval is for use with functions:-argument)

       try
         if type( i, boolean ) and not evalb( evalf( eval( subs( tsubstitutions, i ) ) ) )
         or not is( subs( tsubstitutions, i ) ) then
           WARNING( "constraint `%1' is not satisfied", i );
           WARNING( "subsequent results are not guaranteed" );
         else
           USERINFO( 4, ':-CFSF', `constraint`, i, `is satisfied` );
         end if;
       catch:
         WARNING( "failure in determining whether constraint `%1' is satisfied or not", i );
         WARNING( "subsequent results are not guaranteed" );
       end try;

     end do;

   end if;

   # return the substitutions

   return eval(tsubs);

  end proc;

  # the following procedure constructs the set of substitutions
  # out of the table that processargs returns, and if furthermore
  # the type allows a modification argument, it is also returned
  # in second plase if there is one, otherwise it returns 0;
  # this procedure is used in the evaluation procedures from the
  # submodules (i.e., nth...)

  constructsubs := proc( tformula::table, tsubs::table )

    local tsubstitutions, tmodification;

    # setup the substitutions from the tsubs table

    tsubstitutions := {};

    if assigned( tsubs['parameters'] ) then
      tsubstitutions := tsubstitutions union tsubs['parameters'];
    end if;

    if assigned( tsubs['variable'] ) then
      tsubstitutions := tsubstitutions union { tformula['variable'] = tsubs['variable'] };
    end if;

    # only check for the modification if the type allows it

    if evalb( 'modification' in tformula['type']:-argnames ) then

      if assigned( tsubs['modification'] ) then
        if assigned( tformula['modification'] ) then
          WARNING( "the formula has value `%1' for argument `modification'", eval( tformula['modification'] ) );
          WARNING( "replacing with new value `%1' for this evaluation", eval( tsubs['modification'] ) );
        end if;
        tmodification := tsubs['modification'];
      elif assigned( tformula['modification'] ) then
        tmodification := tformula['modification'];
      else
        tmodification := 0;
      end if;

      # if the modification is of the form tailestimate() or tailestimate('improved'),
      # a modification is computed by calling the internal tailestimate function

      if type( tmodification, function ) and op(0,tmodification) = tailestimate then
        tmodification := tformula['type']:-tailestimate( tformula, op( tmodification ), tsubstitutions );
      end if;

      return ( tsubstitutions, tmodification );

    else
      return ( tsubstitutions );
    end if;

  end proc;

  # the dispatch procedure is used to call the procedure 'fcall' from the
  # submodule specified by the type of the formula 'aformula', when it exists,
  # otherwise an error is returned;
  # when 4 or more arguments are given, we check to see whether args[4] is of
  # type nonnegint; if not, we simply return the unevaluated call using the
  # name specified by 'pname', so the evaluation can be delayed until a true
  # value is known...

  dispatch := proc( pname, fcall, aformula )

    local func, tformula, tsubs;

    tformula := eval(aformula);

    if member( fcall, tformula['type'], 'func' ) then
      if nargs >= 4 and type( eval(args[1]), procedure ) then
        if type( args[4], nonnegint ) then
          tsubs := processargs( pname, fcall, tformula, args[5..-1] );
          func( tformula, args[4], tsubs );
        # elif type( args[4], symbol ) or type( args[4], function ) then
        # elif not type( args[4], complex(extended_numeric) ) then
        else
          if type( aformula, table ) then
            return 'pname'( ':-CFSF:-formula'(aformula), args[4..-1] );
          else
            return 'pname'( aformula, args[4..-1] );
          end if;
        # else
        #   ERROR( "`%1' expects its 2nd argument to be of type `nonnegint', but got `%2'", fcall, args[4] );
        end if;
      else
        func( aformula, args[3..nargs] );
      end if;
    else
      ERROR( "function `%1' is not defined for CFSF formulas of type `%2'", fcall, tformula['type'] );
    end if;

  end proc;

##### ----- RETRIEVING PARTS OF FORMULAS -----

  nthpelement := proc( aformula::formula, n )
    dispatch( procname, ':-nthpelement', aformula, n );
  end proc;

  nthpnumer := proc( aformula::formula, n )
    dispatch( procname, ':-nthpnumer', aformula, n );
  end proc;

  nthpdenom := proc( aformula::formula, n )
    dispatch( procname, ':-nthpdenom', aformula, n );
  end proc;

##### ----- EVALUATING FORMULAS -----

  nthnumer := proc( aformula::formula, n )
    dispatch( procname, ':-nthnumer', aformula, n, args[3..-1] );
  end proc;

  nthdenom := proc( aformula::formula, n )
    dispatch( procname, ':-nthdenom', aformula, n, args[3..-1] );
  end proc;

  nthapprox := proc( aformula::formula, n )
    dispatch( procname, ':-nthapprox', aformula, n, args[3..-1] );
  end proc;

##### ----- WORKING WITH TAILS -----

  nthtail := proc( aformula::formula, n )
    dispatch( procname, ':-nthtail', aformula, n );
  end proc;

  tailestimate := proc( )

    local aformula, tformula;

    # if no formula is specified, just return the unevaluated function
    # call; this way, we can use tailestimate as an expression for the
    # modification part (i.e. delaying its evaluation)

    if nargs = 0 then
      return ':-CFSF:-tailestimate'( );
    elif nargs = 1 and args[1] = 'improved' then
      return ':-CFSF:-tailestimate'( 'improved' );
    end if;

    aformula := args[1];

    if not type( aformula, formula ) then
      ERROR( "`tailestimate' expects its first argument to be of type `formula', but got `%1'", args[1] );
    end if;

    tformula := eval(aformula);

    if member( ':-tailestimate', tformula['type'] ) then
      tformula['type']:-tailestimate( args );
    else
      ERROR( "`tailestimate' is not defined for formulas of type `%1'", tformula['type'] );
    end if;

  end proc;

##### ----- TRANSFORMATIONS AND OTHER FUNCTIONALITY -----

  transform := proc( aformula::formula, form::name )

    local tformula;

    tformula := eval(aformula);

    if member( ':-transform', tformula['type'] ) then
      tformula['type']:-transform( aformula, form );
    else
      ERROR( "no transformations are defined for formulas of type `%1'", tformula['type'] );
    end if;

  end proc;

  # the following function can be used to do a variable
  # substitution, replacing for instance z by z^2 in
  # the definition of the function; the variable itself
  # remains the same

  varsubs := proc( aformula::formula )

    local tformula, tvariable, tsubs;

    if not nargs = 2 then
      ERROR( "`varsubs' expects 2 arguments, but got %1", nargs );
    end if;

    tformula := eval(aformula);
    tvariable := tformula['variable'];
    tsubs := args[2];

    if type( tsubs, equation )
    and evalb( op(1,tsubs) = tvariable )
    and type( op(2,tsubs), exprtype )
    then
      tformula := subs( tvariable = op(2,tsubs), eval(tformula) );
      if type( eval(op(2,tsubs)), symbol ) then
        USERINFO( 4, ':-CFSF', `variable`, tvariable, `changed to`, op(2,tsubs) );
      # elif type( op(2,tsubs), complex(numeric) ) then
      # ERROR( "invalid variable substitution" );
      else
        tformula['variable'] := tvariable;
        USERINFO( 4, ':-CFSF', `occurence of variable`, tvariable, `replaced by`, op(2,tsubs) );
      end if;
      return ':-CFSF:-formula'( tformula );
    else
      ERROR( "`varsubs' expects is second argument to be of the form <variable> = <expression>, but got `%1'", tsubs );
    end if;

  end proc;

  # the equivalent function can be used to check whether two formulas
  # are equivalent to each other, that is, whether they have the same
  # approximants

  equivalent := proc( af1::formula, af2::formula )

    local ti, tm, tcheck, tpos;

    # it suffices to compare the approximants, starting with the begin
    # elements and proceeding with all combinations of general elements
    # (which is their least common multiple)

    tm := lcm( nops( af1['general'] ), nops( af2['general'] ) )
          + max( `if`( assigned( af1['begin'] ), nops( af1['begin'] ), 0 ),
                 `if`( assigned( af2['begin'] ), nops( af2['begin'] ), 0 ) );

    tcheck := map( (evalb@simplify), [ seq( eval( nthapprox( af1, ti ) ) - eval( nthapprox( af2, ti ) ) = 0, ti=0..tm ) ] );

    if member( FAIL, tcheck ) then
      return FAIL;
    elif nops( convert( tcheck, set ) ) = 1 and op( 1, tcheck ) = true then
      return true;
    else
      member( false, tcheck, 'tpos' );
      USERINFO( 4, ':-CFSF', `the approximants first differ at position`, tpoos );
      return false;
    end if;

  end proc;

##### ----- DEFAULT SETTINGS -----

$include "defaults.mpl"

  default := proc( adefault )

    description "function to set default values";

    local items, value;

    if nargs = 0 then
      value := copy(eval(defaults));
      return eval(value);
    end if;

    if nargs > 1 then
      ERROR( "`default' expects at most 1 argument, but got %1", nargs );
    end if;

    items := map( op, {indices(defaults)} );

    if not member( op(1,adefault), items ) then
      ERROR( "`%1' is not a valid default value name", op(1,adefault) );
    elif nops( adefault ) = 1 then
      value := eval(defaults[adefault]);
    elif nops(adefault) > 2 or nops(adefault) = 2 and not type( adefault, equation ) then
      ERROR( "`default' expects its argument to be of the form `<default value name> = <value>'", adefault );
    elif type( defaults[op(1,adefault)], boolean ) and not type( op(2,adefault), boolean ) then
      ERROR( "the default value `%1' must be of type `boolean'", op(1,adefault) );
    elif type( defaults[op(1,adefault)], integer ) and not type( op(2,adefault), integer ) then
      ERROR( "the default value `%1' must be of type `integer'", op(1,adefault) );
    else
      defaults[op(1,adefault)] := op(2,adefault):
      value := defaults[op(1,adefault)];
    end if;

    return eval(value);

  end proc;

##### ----- VERSION INFORMATION -----

  bname::string := "CFSF";
  bdate::string := StringTools:-FormatTime( "%b %d %Y (%T)" ); # we don't include (%T)
  bplatform := kernelopts(version);
  bversion := "";
  btarget := "";

  build := proc()

    local info;

    if nargs = 0 then
      info := bname, " package, Build Date ", bdate;
      if type( bversion, posint ) or type( bversion, string ) and bversion <> "" then
        info := info, ", Build ID ", bversion;
      end if;
      if type ( btarget, string ) and ( btarget <> "" ) then
        info := info, ", Target ", btarget;
      end if;
    elif nargs = 1 then
      if args[1] = 'name' then
        info := bname; #, " package";
      elif args[1] = 'date' then
        info := bdate;
      elif args[1] = 'platform' then
        info := bplatform;
      elif args[1] = 'version' then
        if ( type( bversion, posint ) or type( bversion, string ) and bversion <> "" ) then
          info := bversion;
        else
          WARNING( "no version information available" );
          return;
        end if;
      elif args[1] = 'maketarget' then
        info := btarget;
      else
        ERROR( "`%1' is not a valid argument for `build'", args[1] );
      end if;
    else
      ERROR( "`build' expects 1 argument, but got %1", nargs );
    end if;

    # kernelopts(version) returns a symbol, so we do the same here

    info := convert( cat( info ), symbol );

  end proc;

  # the build information is shown when with(CFSF) is executed
  # (by exporting the user level initialization routine 'init')

  init := proc()
    print( build() );
  end proc;

end module:

##### ===== END OF CFSF MAIN MODULE =====

# the following commands are used to add the CFSF module to
# a repository with name CFSF.lib in the current directory;
# the repository is also optimized and garbage collection
# is performed to remove unused entries
# (see also ?repository,management)

try
  savelibname := cat( currentdir(), "/CFSF.lib" ):
  try march( 'create', savelibname, 500 ) catch: end try:
  savelib( 'CFSF' ):
  map( march, ['gc','reindex','pack'], savelibname ):
end try:
