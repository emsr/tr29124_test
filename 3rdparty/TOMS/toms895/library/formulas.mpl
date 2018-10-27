##### ===================================
##### ===== TABLE OF KNOWN FORMULAS =====
##### ===================================

# here, we set the variable 'defaults[add_formulas]' to 'true',
# so that we can use the 'create' function to add formulas to
# the formulas table

defaults['add_formulas'] := true;

# for each function or constant, a file exists which contains
# 'create'-statements (including a 'label' argument) for all
# formulas that are mentioned in the book; these files are included
# below, so all these formulas will be added to the formulas table

# a warning message will be issued if a formula is stored in an object
# that is not known as a local object to the main CFSF module (which
# is the case when its name does not occur in the 'local' or 'export'
# list of the main CFSF module); this warning can be suppressed by
# setting 'warnlevel' to zero, but remember that this object will
# not be accessible by its object name

# formulas for chapter CN

$include "../functions/CN/archimedes/archimedes.mpl"
$include "../functions/CN/eulersnumber/eulersnumber.mpl"
$include "../functions/CN/powerandroot/powerandroot.mpl"
$include "../functions/CN/naturallogarithm/naturallogarithm.mpl"
$include "../functions/CN/pythagoras/pythagoras.mpl"
$include "../functions/CN/eulersconstant/eulersconstant.mpl"
$include "../functions/CN/goldenratio/goldenratio.mpl"
$include "../functions/CN/rabbit/rabbit.mpl"
$include "../functions/CN/apery/apery.mpl"
$include "../functions/CN/catalan/catalan.mpl"
$include "../functions/CN/gompertz/gompertz.mpl"

# formulas for chapter EF

$include "../functions/EF/exp/exp.mpl"
$include "../functions/EF/ln/ln.mpl"

$include "../functions/EF/sin/sin.mpl"
$include "../functions/EF/cos/cos.mpl"
$include "../functions/EF/tan/tan.mpl"

$include "../functions/EF/arcsin/arcsin.mpl"
$include "../functions/EF/arccos/arccos.mpl"
$include "../functions/EF/arctan/arctan.mpl"

$include "../functions/EF/sinh/sinh.mpl"
$include "../functions/EF/cosh/cosh.mpl"
$include "../functions/EF/tanh/tanh.mpl"
$include "../functions/EF/coth/coth.mpl"

$include "../functions/EF/arcsinh/arcsinh.mpl"
$include "../functions/EF/arccosh/arccosh.mpl"
$include "../functions/EF/arctanh/arctanh.mpl"

$include "../functions/EF/pow/pow.mpl"

# formulas for Chapter GA

$include "../functions/GA/binet/binet.mpl"
$include "../functions/GA/polygamma/polygamma.mpl"
$include "../functions/GA/trigamma/trigamma.mpl"
$include "../functions/GA/tetragamma/tetragamma.mpl"
$include "../functions/GA/incompletegamma/incompletegamma.mpl"

# formulas for Chapter ER

$include "../functions/ER/comperror/comperror.mpl"
$include "../functions/ER/error/error.mpl"
$include "../functions/ER/fresnel/fresnel.mpl"
$include "../functions/ER/repint/repint.mpl"

# formulas for Chapter EX

$include "../functions/EX/expintegrals/expintegrals.mpl"
$include "../functions/EX/related/related.mpl"

# formulas for Chapter HY

$include "../functions/HY/hypergeometric.mpl"

# formulas for Chapter CH

$include "../functions/CH/kummer/kummer.mpl"
$include "../functions/CH/confluent/confluent.mpl"
$include "../functions/CH/confluentlimit/confluentlimit.mpl"
$include "../functions/CH/parabolic/parabolic.mpl"
$include "../functions/CH/whittaker/whittaker.mpl"

# formulas for Chapter BS

$include "../functions/BS/bessel/bessel.mpl"
$include "../functions/BS/modbessel/modbessel.mpl"

# formulas for Chapter SM

$include "../functions/SM/normal/normal.mpl"
$include "../functions/SM/repeated/repeated.mpl"
$include "../functions/SM/gamma_chisquare/gamma_chisquare.mpl"
$include "../functions/SM/beta_f_t/beta_f_t.mpl"

# formulas for Chapter QH

$include "../functions/QH/qhyper.mpl"

# NOTE : if you want another initial formulas table,
# you can define it more explicitly too; e.g.
#
# formula := table([
#   "11.1.1" = EF.exp.power.01,
#   "11.1.2" = EF.exp.sfrac.01,
#   "11.1.3" = EF.exp.cfrac.01,
#   "11.1.4" = EF.exp.tfrac.01
# ]);

# finally, we reset the variable 'defaults[add_formulas]' to 'false',
# its original value

defaults['add_formulas'] := false;

USERINFO( 3, ':-CFSF', nops(op(op(formulas))), `formulas are available` );
