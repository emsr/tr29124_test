/* quadrature/qaws_integration_table.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
 * Copyright (C) 2016 Edward Smith-Rowland
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#ifndef QAWS_INTEGRATION_TABLE_H
#define QAWS_INTEGRATION_TABLE_H 1

#include <array>

namespace __gnu_test
{

template<typename _Tp>
  struct qaws_integration_table
  {
    _Tp alpha;
    _Tp beta;
    int mu;
    int nu;
    std::array<_Tp, 25> ri;
    std::array<_Tp, 25> rj;
    std::array<_Tp, 25> rg;
    std::array<_Tp, 25> rh;

    qaws_integration_table(_Tp alpha_in, _Tp beta_in, int mu_in, int nu_in);
    void set(_Tp alpha_in, _Tp beta_in, int mu_in, int nu_in);
    void initialise();
  };

template<typename _Tp>
  qaws_integration_table<_Tp>::qaws_integration_table(_Tp alpha_in, _Tp beta_in, int mu_in, int nu_in)
  : alpha(alpha_in),
    beta(beta_in),
    mu(mu_in),
    nu(nu_in)
  {
    if (this->alpha < _Tp{-1})
      std::__throw_domain_error ("alpha must be greater than -1.0");
    if (this->beta < _Tp{-1})
      std::__throw_domain_error ("beta must be greater than -1.0");
    if (this->mu != 0 && this->mu != 1)
      std::__throw_domain_error ("mu must be 0 or 1");
    if (this->nu != 0 && this->nu != 1)
      std::__throw_domain_error ("nu must be 0 or 1");

    this->initialise();
  }


template<typename _Tp>
  void
  qaws_integration_table<_Tp>::set(_Tp alpha_in, _Tp beta_in, int mu_in, int nu_in)
  {
    if (alpha_in < _Tp{-1})
      std::__throw_domain_error ("alpha must be greater than -1.0");
    if (beta_in < _Tp{-1})
      std::__throw_domain_error ("beta must be greater than -1.0");
    if (mu_in != 0 && mu_in != 1)
      std::__throw_domain_error ("mu must be 0 or 1");
    if (nu_in != 0 && nu_in != 1)
      std::__throw_domain_error ("nu must be 0 or 1");

    this->alpha = alpha_in;
    this->beta = beta_in;
    this->mu = mu_in;
    this->nu = nu_in;

    this->initialise();
  }

template<typename _Tp>
  void
  qaws_integration_table<_Tp>::initialise()
  {
    const auto alpha_p1 = this->alpha + _Tp{1};
    const auto beta_p1 = this->beta + _Tp{1};

    const auto alpha_p2 = this->alpha + _Tp{2};
    const auto beta_p2 = this->beta + _Tp{2};

    const auto r_alpha = std::pow(_Tp{2}, alpha_p1);
    const auto r_beta = std::pow(_Tp{2}, beta_p1);

    this->ri[0] = r_alpha / alpha_p1;
    this->ri[1] = this->ri[0] * this->alpha / alpha_p2;
    auto an = _Tp{2};
    auto anm1 = _Tp{1};
    for (size_t i = 2; i < this->ri.size(); ++i)
      {
	this->ri[i] = -(r_alpha + an * (an - alpha_p2) * this->ri[i - 1])
          / (anm1 * (an + alpha_p1));
	anm1 = an;
	an += _Tp{1};
      }

    this->rj[0] = r_beta / beta_p1;
    this->rj[1] = this->rj[0] * this->beta / beta_p2;
    an = _Tp{2};
    anm1 = _Tp{1};
    for (size_t i = 2; i < this->rj.size(); ++i)
      {
	this->rj[i] = -(r_beta + an * (an - beta_p2) * this->rj[i - 1])
          / (anm1 * (an + beta_p1));
	anm1 = an;
	an += _Tp{1};
      }
    for (size_t i = 1; i < this->rj.size(); i += 2)
      this->rj[i] *= _Tp{-1};

    this->rg[0] = -this->ri[0] / alpha_p1;
    this->rg[1] = -this->rg[0] - 2.0 * r_alpha / (alpha_p2 * alpha_p2);
    an = _Tp{2};
    anm1 = _Tp{1};
    for (size_t i = 2; i < this->rg.size(); ++i)
      {
	this->rg[i] = -(an * (an - alpha_p2) * this->rg[i - 1] - an * this->ri[i - 1]
                  + anm1 * this->ri[i]) / (anm1 * (an + alpha_p1));
	anm1 = an;
	an += _Tp{1};
      }

    this->rh[0] = -this->rj[0] / beta_p1;
    this->rh[1] = -this->rh[0] - _Tp{2} * r_beta / (beta_p2 * beta_p2);
    an = _Tp{2};
    anm1 = _Tp{1};
    for (size_t i = 2; i < this->rh.size(); ++i)
      {
	this->rh[i] = -(an * (an - beta_p2) * this->rh[i - 1] - an * this->rj[i - 1]
                  + anm1 * this->rj[i]) / (anm1 * (an + beta_p1));
	anm1 = an;
	an += _Tp{1};
      }
    for (size_t i = 1; i < this->rh.size(); i += 2)
      this->rh[i] *= _Tp{-1};
  }

} // namespace

#endif // QAWS_INTEGRATION_TABLE_H
