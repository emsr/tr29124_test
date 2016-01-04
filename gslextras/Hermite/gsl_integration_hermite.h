// gsl_integration_hermite.h

//*----------------------------------------------------------------------*
//* This program is free software: you can redistribute it and/or modify *
//* it under the terms of the GNU General Public License as published by *
//* the Free Software Foundation, either version 3 of the License, or    *
//* (at your option) any later version.                                  *
//*                                                                      *
//* This program is distributed in the hope that it will be useful,      *
//* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
//* GNU General Public License for more details.                         *
//*                                                                      *
//* You should have received a copy of the GNU General Public License    *
//* along with this program. If not, see <http://www.gnu.org/licenses/>. *
//*----------------------------------------------------------------------*

//*----------------------------------------------------------------------*
//* Copyright 2014 Konrad Griessinger                                    *
//* (konradg(at)gmx.net)                                                 *
//*----------------------------------------------------------------------*



#ifndef __GSL_INTEGRATION_HERMITE_H__
#define __GSL_INTEGRATION_HERMITE_H__

#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

double gsl_integration_hermite_prob(const int n, const gsl_function *func);
double gsl_integration_hermite_phys(const int n, const gsl_function *func);

__END_DECLS

#endif /* __GSL_INTEGRATION_HERMITE_H__ */
