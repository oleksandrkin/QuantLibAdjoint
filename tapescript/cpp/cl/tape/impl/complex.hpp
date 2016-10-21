/*
Copyright (C) 2015-present CompatibL

Performance test results and finance-specific examples are available at:

http://www.tapescript.org

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef cl_tape_impl_complex_hpp
#define cl_tape_impl_complex_hpp

#include <cl/tape/impl/tape_fwd.hpp>
#include <complex>

namespace ext = CppAD;

namespace std
{
    /// <summary>The complex type based on tape double.</summary>
    template<typename Inner>
    class complex<cl::tape_wrapper<Inner>>
    {
    public:

        // The type holder type
        typedef cl::tape_wrapper<Inner> tape_type;

        // parse ad base
        typedef typename
            cl::remove_ad<typename tape_type::value_type>::type ad_base_type;

        typedef tape_type value_type;
        typedef tape_type real_type;
        typedef complex<tape_type> complex_type;
        typedef complex<ad_base_type> complex_double;
        typedef cl::tape_inner_type<complex_double > complex_based_type;
        typedef std::complex<cl::tape_inner_type<ad_base_type> > real_based_type;

        enum Complex_Mode { None = 0, RealBase = (1 << 1), ComplBase = (1 << 2) };

        static const Complex_Mode default_mode = ComplBase;

        // If type is initialized by certain values this is real base type
        complex(real_type const& real, real_type const& imag = 0.0, Complex_Mode mode = default_mode)
            : real_base_()
            , complex_base_()
            , mode_(mode)
        {
            if (mode_ == RealBase)
            {
                real_base_ = real_based_type(real.value(), imag.value());
            }
            else
            {
                bool is_variable = ext::Variable(cl::tapescript::cvalue(real))
                    || ext::Variable(cl::tapescript::cvalue(imag));

                if (is_variable)
                {
                    mode_ = RealBase;
                    real_base_ = real_based_type(real.value(), imag.value());
                }
                else
                {
                    complex_base_ = complex_double(ext::Value(real.value()), ext::Value(imag.value()));
                }
            }
        }
        
        complex(double const& real, double const& imag = 0.0, Complex_Mode mode = default_mode)
            : real_base_()
            , complex_base_()
            , mode_(mode)
        {
            if (mode_ == RealBase)
            {
                real_base_ = real_based_type((tape_type)real, (tape_type)imag);
            }
            else
            {
                bool is_variable = ext::Variable(cl::tapescript::cvalue((tape_type)real))
                    || ext::Variable(cl::tapescript::cvalue((tape_type)imag));

                if (is_variable)
                {
                    mode_ = RealBase;
                    real_base_ = real_based_type((tape_type)real, (tape_type)imag);
                }
                else
                {
                    complex_base_ = complex_double(ext::Value(((tape_type)real).value()), ext::Value(((tape_type)imag).value()));
                }
            }
        }

        // This will be called if vector is resized
        complex()
            : real_base_()
            , complex_base_()
            , mode_(default_mode)
        {    }

        complex(complex const& other)
            : real_base_(other.real_base_)
            , complex_base_(other.complex_base_)
            , mode_(other.mode_)
        {    }

        template<class Ty>
        explicit complex(complex<Ty> const & other, Complex_Mode mode = default_mode)
            : complex(other.real(), other.imag(), mode)
        {     }

        explicit complex(const complex_based_type & other)
            : real_base_()
            , complex_base_(other)
            , mode_(ComplBase)
        {    }

        inline void check_mode(Complex_Mode mode) const
        {
            if (mode != mode_)
                cl::throw_("Wrong tape mode.");
        }

        inline void set_mode(Complex_Mode mode)
        {
            if (mode == mode_)
                return;

            if (mode == RealBase)
            {
                real_base_ = to_real_base();
                mode_ = RealBase;
            }
            else
            {
                complex_base_ = to_complex_base();
                mode_ = ComplBase;
            }
        }

        inline real_type real(real_type const& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_.real(right.value());
            }
            else
            {
                complex_base_ = complex_double(ext::Value(right.value()), ext::Value(complex_base_).imag());
            }
            return right;
        }

        inline real_type imag(const real_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_.imag(right.value());
            }
            else
            {
                complex_base_ = complex_double(ext::Value(complex_base_).real(), ext::Value(right.value()));
            }
            return right;
        }

        inline real_type real() const
        {
            if (mode_ == RealBase)
            {
                return real_base_.real();
            }
            return ext::Value(complex_base_).real();
        }

        inline real_type imag() const
        {
            if (mode_ == RealBase)
            {
                return real_base_.imag();
            }
            return ext::Value(complex_base_).imag();
        }

        inline complex_type& operator=(complex_type const& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_ = right.to_real_base();
            }
            else
            {
                complex_base_ = right.to_complex_base();
            }

            return *this;
        }

        inline complex_type& operator=(double const& right)
        {
            if (mode_ == RealBase)
            {
                real_base_ = right;
            }
            if (mode_ == ComplBase)
            {
                complex_base_ = right;
            }
            return (*this);
        }

        inline complex_type& operator=(const real_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_ = right.value();
            }
            else
            {
                complex_base_ = ext::Value(right.value());
            }
            return (*this);
        }

        inline complex_type& operator+=(const real_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_ += right.value();
            }
            else
            {
                complex_base_ += ext::Value(right.value());
            }
            return (*this);
        }

        inline complex_type& operator-=(const real_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_ -= right.value();
            }
            else
            {
                complex_base_ -= ext::Value(right.value());
            }
            return (*this);
        }

        inline complex_type& operator*=(const real_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_ *= right.value();
            }
            else
            {
                complex_base_ *= ext::Value(right.value());
            }
            return (*this);
        }

        inline complex_type& operator/=(const real_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_ /= right.value();
            }
            else
            {
                complex_base_ /= ext::Value(right.value());
            }
            return (*this);
        }

        inline complex_type& operator+=(const complex_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_ += right.to_real_base();
            }
            else
            {
                complex_base_ += right.to_complex_base();
            }
            return (*this);
        }

        inline complex_type& operator-=(const complex_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_ -= right.to_real_base();
            }
            else
            {
                complex_base_ -= right.to_complex_base();
            }
            return (*this);
        }

        inline complex_type& operator*=(const complex_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
                real_base_ *= right.to_real_base();
            }
            else
            {
                complex_base_ *= right.to_complex_base();
            }
            return (*this);
        }

        inline complex_type& operator/=(const complex_type& right)
        {
            if (adjust_mode(right) == RealBase)
            {
#define CL_OPERATOR_DIV_EQ_FIXED_NOT
#if defined CL_OPERATOR_DIV_EQ_FIXED
                real_base_ /= right.to_real_base();
#else
                real_based_type right_base = right.to_real_base();
                cl::tape_inner_type<ad_base_type> x1 = real_base_.real()
                    , y1 = real_base_.imag()
                    , x2 = right_base.real()
                    , y2 = right_base.imag()
                    , norm2 = pow(x2, 2) + pow(y2, 2);
                real_base_ = real_based_type(
                    (x1 * x2 + y1 * y2) / norm2
                    , (y1 * x2 - x1 * y2) / norm2);
                // cl::throw_("Can't use operator: " __FUNCSIG__);
#endif
            }
            else
            {
                complex_base_ /= right.to_complex_base();
            }
            return (*this);
        }

        inline complex_type& __rep()
        {
            return (*this);
        }

        // Adjust the mode that should be used for calculation.
        Complex_Mode adjust_mode(real_type const& right)
        {
            if (mode_ == RealBase)
                return mode_;

            if (ext::Variable(right.value()))
            {
                cl::CheckParameter(complex_base_);
                real_base_ = real_based_type(ext::Value(complex_base_).real(), ext::Value(complex_base_).imag());
                mode_ = RealBase;
            }
            return mode_;
        }

        // Adjust the mode that should be used for calculation.
        Complex_Mode adjust_mode(complex_type const& right)
        {
            if (mode_ == right.mode_)
                return mode_;

            if (mode_ == RealBase)
            {
                if (ext::Variable(right.complex_base_))
                {
                    auto real_temp = real_base_.real();
                    cl::CheckParameter(real_temp);
                    auto imag_temp = real_base_.imag();
                    cl::CheckParameter(imag_temp);
                    complex_base_ = complex_double(ext::Value(real_base_.real()), ext::Value(real_base_.imag()));
                    mode_ = ComplBase;
                }
                return mode_;
            }

            if (ext::Variable(right.real_base_.real()) || ext::Variable(right.real_base_.imag()))
            {
                cl::CheckParameter(complex_base_);
                real_base_ = real_based_type(ext::Value(complex_base_).real(), ext::Value(complex_base_).imag());
                mode_ = RealBase;
            }
            return mode_;
        }

        complex_based_type& complex_base()
        {
            set_mode(ComplBase);
            return complex_base_;
        }

        const complex_based_type& complex_base() const
        {
            check_mode(ComplBase);
            return complex_base_;
        }

        complex_based_type to_complex_base() const
        {
            if (mode_ == RealBase)
            {
                return complex_double(ext::Value(real_base_.real()), ext::Value(real_base_.imag()));
            }
            return complex_base_;
        }

        real_based_type to_real_base() const
        {
            if (mode_ == RealBase)
            {
                return real_base_;
            }
            return real_based_type(real(), imag());
        }

        inline complex<tape_type> operator - () const
        {
            return complex<tape_type>(0.0) - *this;
        }

        // Used when tape is recording for real or image part of complex value.
        real_based_type real_base_;
        // Used when tape is recording for the complex value as one variable.
        complex_based_type complex_base_;
        Complex_Mode mode_;
    };
}

namespace cl
{
    // returns true iff the absolute value of lhs is greater than or equal rhs
    template <typename Inner>
    inline bool abs_geq(const std::complex<tape_wrapper<Inner>>& lhs, const double& rhs)
    {
        typedef typename
            cl::remove_ad<typename tape_wrapper<Inner>::value_type>::type inner_base_type;
        if (lhs.mode_
            == std::complex<tape_wrapper<Inner>>::RealBase)
        {
            return std::abs(lhs) >= rhs;
        }
        else
        {
            return ext::abs_geq(lhs.complex_base_, cl::tape_inner_type<std::complex<inner_base_type>>(rhs));
        }
    }
}

namespace std
{
    // Arithmetics for std::complex<cl::tape_double> and cl::tape_double

    inline std::complex<cl::tape_double> operator+(
        const std::complex<cl::tape_double>& lhs, cl::tape_double rhs)
    {
        complex<cl::tape_double> temp = lhs;
        return temp += rhs;
    }

    inline std::complex<cl::tape_double> operator-(
        const std::complex<cl::tape_double>& lhs, cl::tape_double rhs)
    {
        complex<cl::tape_double> temp = lhs;
        return temp -= rhs;
    }

    inline std::complex<cl::tape_double> operator*(
        const std::complex<cl::tape_double>& lhs, cl::tape_double rhs)
    {
        complex<cl::tape_double> temp = lhs;
        return temp *= rhs;
    }

    inline std::complex<cl::tape_double> operator/(
        const std::complex<cl::tape_double>& lhs, cl::tape_double rhs)
    {
        complex<cl::tape_double> temp = lhs;
        return temp /= rhs;
    }

    inline std::complex<cl::tape_double> operator+(
        cl::tape_double lhs, const std::complex<cl::tape_double>& rhs)
    {
        complex<cl::tape_double> temp = rhs;
        return temp += lhs;
    }

    inline std::complex<cl::tape_double> operator-(
        cl::tape_double lhs, const std::complex<cl::tape_double>& rhs)
    {
        complex<cl::tape_double> temp = lhs;
        return temp -= rhs;
    }

    // Arithmetics for std::complex<cl::tape_double> and double
    inline std::complex<cl::tape_double>
    operator*(cl::tape_double lhs, const std::complex<cl::tape_double>& rhs)
    {
        complex<cl::tape_double> temp = rhs;
        return temp *= lhs;
    }

    inline std::complex<cl::tape_double>
    operator/(cl::tape_double lhs, const std::complex<cl::tape_double>& rhs)
    {
#if defined CL_TAPE_COMPLEX_ENABLED
        complex<cl::tape_double> temp = lhs;
        return temp /= rhs;
#else
        return std::complex<cl::tape_double>(lhs, 0.0) / rhs;
#endif
    }

    // Arithmetics for std::complex<cl::tape_double> and double
#if defined CL_TAPE_COMPLEX_ENABLED
    complex<cl::tape_double> inline pow_(
        complex<cl::tape_double> const &_Left
        , complex<cl::tape_double> const &_First
        , int _Right, bool _Even)
    {
        if (_Right <= 1)
            return _Even ? _Left : _Left*_First;

        return pow_(_Left * _Left, _First, _Right / 2, _Even);
    }

    //template <typename Right>
    inline complex<cl::tape_double>
    pow(complex<cl::tape_double> const & _Left, int _Right)
    {
        if (_Right == 0)
            return complex<cl::tape_double>(1.0);

        return _Right == 1 ? _Left : pow_(_Left * _Left, _Left
            , _Right / 2, (_Right / 2) * 2 == _Right);
    }

#endif

    inline std::complex<cl::tape_double>
    operator+(const std::complex<cl::tape_double>& lhs, double rhs)
    { return lhs + cl::tape_double(rhs); }

    inline std::complex<cl::tape_double>
    operator-(const std::complex<cl::tape_double>& lhs, double rhs)
    { return lhs - cl::tape_double(rhs); }

    inline std::complex<cl::tape_double>
    operator*(const std::complex<cl::tape_double>& lhs, double rhs)
    { return lhs * cl::tape_double(rhs); }

    inline std::complex<cl::tape_double>
    operator/(const std::complex<cl::tape_double>& lhs, double rhs)
    { return lhs / cl::tape_double(rhs); }

    inline std::complex<cl::tape_double>
    operator+(double lhs, const std::complex<cl::tape_double>& rhs)
    { return cl::tape_double(lhs) + rhs; }

    inline std::complex<cl::tape_double>
    operator-(double lhs, const std::complex<cl::tape_double>& rhs)
    { return cl::tape_double(lhs) - rhs; }

    inline std::complex<cl::tape_double>
    operator*(double lhs, const std::complex<cl::tape_double>& rhs)
    { return cl::tape_double(lhs) * rhs; }

    inline std::complex<cl::tape_double>
    operator/(double lhs, const std::complex<cl::tape_double>& rhs)
    { return cl::tape_double(lhs) / rhs; }

#if defined CL_TAPE_COMPLEX_ENABLED

#define CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(FUN)                                        \
    template<typename Inner>                                                                            \
    inline std::complex<cl::tape_wrapper<Inner>> FUN(const std::complex<cl::tape_wrapper<Inner>>& lhs)    \
    {                                                                                   \
        if (lhs.mode_ == std::complex<cl::tape_wrapper<Inner>>::RealBase)                                \
        {                                                                               \
            return std::FUN<cl::tape_wrapper<Inner>>(lhs);                                               \
        }                                                                               \
        return std::complex<cl::tape_wrapper<Inner>>(FUN(lhs.complex_base_));                            \
    }

    CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(sqrt)
    CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(exp)
    CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(log)
    CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(log10)
    CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(sin)
    CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(cos)
    CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(sinh)
    CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(cosh)
    //CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(acos)
    //CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(asin)
    //CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC(atan)

    inline std::complex<cl::tape_double> tan(const std::complex<cl::tape_double>& lhs)
    {
        if (lhs.mode_ == std::complex<cl::tape_double>::RealBase)
        {
            return std::tan<cl::tape_double>(lhs);
        }
        return std::complex<cl::tape_double>(
            sin(lhs.complex_base_) / cos(lhs.complex_base_));
    }

    inline std::complex<cl::tape_double> tanh(const std::complex<cl::tape_double>& lhs)
    {
        if (lhs.mode_ == std::complex<cl::tape_double>::RealBase)
        {
            return std::tanh<cl::tape_double>(lhs);
        }
        return std::complex<cl::tape_double>(
            sinh(lhs.complex_base_) / cosh(lhs.complex_base_));
    }

    inline complex<cl::tape_double> pow(
        complex<cl::tape_double> const &_Left, cl::tape_double _Right)
    {
        return exp(_Right * log(_Left));
    }

    inline complex<cl::tape_double> pow(
        cl::tape_double const &_Left, complex<cl::tape_double> const &_Right)
    {
        return exp(_Right * log(complex<cl::tape_double>(_Left)));
    }

    inline complex<cl::tape_double> pow(
        complex<cl::tape_double> const &_Left
        , complex<cl::tape_double> const &_Right)
    {
        return exp(_Right * log(_Left));
    }


#undef CL_TAPE_DOUBLE_COMPLEX_GENERIC_FUNC

#endif

    //!! Providing implementation causes compilation error due to NaN not being defined for Real, to be resolved
}

#endif // cl_tape_impl_complex_hpp
