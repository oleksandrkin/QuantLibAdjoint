/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

#include <ql/PricingEngines/Vanilla/analyticdividendeuropeanengine.hpp>
#include <ql/PricingEngines/blackformula.hpp>
#include <ql/Processes/blackscholesprocess.hpp>

namespace QuantLib {

    void AnalyticDividendEuropeanEngine::calculate() const {

        QL_REQUIRE(arguments_.exercise->type() == Exercise::European,
                   "not an European option");

        boost::shared_ptr<StrikedTypePayoff> payoff =
            boost::dynamic_pointer_cast<StrikedTypePayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-striked payoff given");

        boost::shared_ptr<BlackScholesProcess> process =
            boost::dynamic_pointer_cast<BlackScholesProcess>(
                                                arguments_.stochasticProcess);
        QL_REQUIRE(process, "Black-Scholes process required");

        Date settlementDate = process->riskFreeRate()->referenceDate();
        Real riskless = 0.0;
        Size i;
        for (i=0; i<arguments_.dividends.size(); i++)
            if (arguments_.dividendDates[i] >= settlementDate)
                riskless += arguments_.dividends[i] *
                    process->riskFreeRate()
                    ->discount(arguments_.dividendDates[i]);
        Real spot = process->stateVariable()->value() - riskless;

        DiscountFactor dividendDiscount =
            process->dividendYield()->discount(
                                            arguments_.exercise->lastDate());
        DiscountFactor riskFreeDiscount =
            process->riskFreeRate()->discount(
                                            arguments_.exercise->lastDate());
        Real forwardPrice = spot * dividendDiscount / riskFreeDiscount;

        Real variance =
            process->blackVolatility()->blackVariance(
                                              arguments_.exercise->lastDate(),
                                              payoff->strike());

        BlackFormula black(forwardPrice, riskFreeDiscount, variance, payoff);

        results_.value = black.value();
        results_.delta = black.delta(spot);
        results_.gamma = black.gamma(spot);

        DayCounter rfdc  = process->riskFreeRate()->dayCounter();
        DayCounter voldc = process->blackVolatility()->dayCounter();
        Time t = voldc.yearFraction(
                                  process->blackVolatility()->referenceDate(),
                                  arguments_.exercise->lastDate());
        results_.vega = black.vega(t);

        Real delta_theta = 0.0, delta_rho = 0.0;
        for (i = 0; i < arguments_.dividends.size(); i++) {
            Date d = arguments_.dividendDates[i];
            if (d >= settlementDate) {
                delta_theta -= arguments_.dividends[i] *
                  process->riskFreeRate()->zeroRate(d,rfdc,Continuous,Annual)*
                  process->riskFreeRate()->discount(d);
                Time t = rfdc.yearFraction(
                    process->riskFreeRate()->referenceDate(), d);
                delta_rho += arguments_.dividends[i] * t *
                             process->riskFreeRate()->discount(t);
            }
        }
        t = rfdc.yearFraction(process->riskFreeRate()->referenceDate(),
                              arguments_.exercise->lastDate());
        try {
            results_.theta = black.theta(spot, t) +
                             delta_theta * black.delta(spot);
        } catch (Error&) {
            results_.theta = Null<Real>();
        }

        results_.rho = black.rho(t) +
                       delta_rho * black.delta(spot);
    }

}

