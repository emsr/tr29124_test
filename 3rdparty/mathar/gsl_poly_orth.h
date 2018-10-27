/* poly/gsl_poly_orth.h
*
* 2006-09-09 Richard J. Mathar, Leiden Observatory
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef __GSL_POLY_ORTH_H__
#define __GSL_POLY_ORTH_H__


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

#include <stdlib.h>

int gsl_poly_dd_init_orthH (double dd[], double xa[], const size_t size) ;
int gsl_poly_dd_init_orthL (double dd[], double xa[], const size_t size) ;
int gsl_poly_dd_init_orthT (double dd[], double xa[], const size_t size) ;

__END_DECLS

#endif /* __GSL_POLY_ORTH_H__ */
