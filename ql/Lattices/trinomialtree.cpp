/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/reference/license.html>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/Lattices/trinomialtree.hpp>

namespace QuantLib {

    TrinomialTree::TrinomialTree(
                           const boost::shared_ptr<StochasticProcess>& process,
                           const TimeGrid& timeGrid,
                           bool isPositive)
    : Tree(timeGrid.size()), dx_(1, 0.0), timeGrid_(timeGrid) {
        x0_ = process->x0();

        Size nTimeSteps = timeGrid.size() - 1;
        Integer jMin = 0;
        Integer jMax = 0;

        for (Size i=0; i<nTimeSteps; i++) {
            Time t = timeGrid[i];
            Time dt = timeGrid.dt(i);

            //Variance must be independent of x
            Real v2 = process->variance(t, 0.0, dt);
            Volatility v = std::sqrt(v2);
            dx_.push_back(v*std::sqrt(3.0));

            boost::shared_ptr<TrinomialBranching> branching(
                                                    new TrinomialBranching());
            for (Integer j=jMin; j<=jMax; j++) {
                Real x = x0_ + j*dx_[i];
                Real m = process->expectation(t, x, dt);
                Integer temp = Integer(std::floor((m-x0_)/dx_[i+1] + 0.5));

                if (isPositive) {
                    while (x0_+(temp-1)*dx_[i+1]<=0) {
                        temp++;
                    }
                }

                branching->k_.push_back(temp);
                Real e = m - (x0_ + temp*dx_[i+1]);
                Real e2 = e*e;
                Real e3 = e*std::sqrt(3.0);

                branching->probs_[0].push_back((1.0 + e2/v2 - e3/v)/6.0);
                branching->probs_[1].push_back((2.0 - e2/v2)/3.0);
                branching->probs_[2].push_back((1.0 + e2/v2 + e3/v)/6.0);
            }
            branchings_.push_back(branching);

            const std::vector<Integer>& k = branching->k_;
            jMin = *std::min_element(k.begin(), k.end()) - 1;
            jMax = *std::max_element(k.begin(), k.end()) + 1;
        }

    }

    Real TrinomialTree::underlying(Size i, Size index) const {
        if (i==0) return x0_;
        const std::vector<Integer>& k = branchings_[i-1]->k_;
        Integer jMin = *std::min_element(k.begin(), k.end()) - 1;
        return x0_ + (jMin*1.0 + index*1.0)*dx(i);
    }

    Size TrinomialTree::size(Size i) const {
        if (i==0) return 1;
        const std::vector<Integer>& k = branchings_[i-1]->k_;
        Integer jMin = *std::min_element(k.begin(), k.end()) - 1;
        Integer jMax = *std::max_element(k.begin(), k.end()) + 1;
        return jMax - jMin + 1;
    }

}

