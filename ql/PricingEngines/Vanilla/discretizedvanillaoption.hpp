/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2004 StatPro Italia srl

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

/*! \file discretizedvanillaoption.hpp
    \brief discretized vanilla option
*/

#ifndef quantlib_discretized_vanilla_option_h
#define quantlib_discretized_vanilla_option_h

#include <ql/discretizedasset.hpp>
#include <ql/Lattices/bsmlattice.hpp>
#include <ql/Pricers/singleassetoption.hpp>
#include <ql/Instruments/vanillaoption.hpp>

namespace QuantLib {

    class DiscretizedVanillaOption : public DiscretizedAsset {
      public:
        DiscretizedVanillaOption(const VanillaOption::arguments& args)
        : arguments_(args) {}

        void reset(Size size);

        std::vector<Time> mandatoryTimes() const {
            return arguments_.stoppingTimes;
        }
      protected:
        void postAdjustValuesImpl();
      private:
        void applySpecificCondition();
        VanillaOption::arguments arguments_;
    };

}





#endif
