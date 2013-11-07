// Copyright (C) 2013-2014 by Paul Cazeaux
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


%define BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(CLASS)
    %template(CLASS ## _float32_float32)
        CLASS<float, float>;
    %template(CLASS ## _float32_complex64)
        CLASS<float, std::complex<float> >;
    %template(CLASS ## _complex64_complex64)
        CLASS<std::complex<float>, std::complex<float> >;

    %template(CLASS ## _float64_float64)
        CLASS<double, double>;
    %template(CLASS ## _float64_complex128)
        CLASS<double, std::complex<double> >;
    %template(CLASS ## _complex128_complex128)
        CLASS<std::complex<double>, std::complex<double> >
%enddef // BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT




%module core_periodic_laplace


// %include "assembly/laplace_3d_periodic_potential_operators.i"
%include "assembly/laplace_3d_periodic_operators.i"