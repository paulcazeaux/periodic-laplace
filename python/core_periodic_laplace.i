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


%module core_periodic_laplace
%{
#define SWIG_FILE_WITH_INIT
#include <dune/common/exceptions.hh>
#include <complex>
%}


%include "numpy.i"

// Import docstring macros
%include "docstrings.i"

// Define commonly used typemaps
%include "armadillo.i"
%include "auto_ptr.i"
%include "exception.i"
%include "bempp_shared_ptr.i"
%include "std_string.i"
%include "std_complex.i"

// Useful macros
%include "macros.i"


// Import macros for explicit template instantiations
%include "template_instantiations_basis_result.i"
%include "template_instantiations_value.i"



// End of auto_ptr typemaps


// Make commonly used typedefs known to Swig
%inline %{
    namespace Bempp
    {
        typedef double ctype;
    }
%}
%include "grid/geometry_type.hpp"

// Wrap Bempp components

// Common
%include "common/scalar_traits.i"
// Grid
%include "grid/geometry.i"


// Assembly
%include "assembly/python_surface_normal_dependent_functor.i"
%include "assembly/python_surface_normal_independent_functor.i"
%include "bempp.swg"

%include "assembly/integrate.i"
%include "assembly/laplace_3d_periodic_operators.i"
%include "assembly/laplace_3d_periodic_potential_operators.i"