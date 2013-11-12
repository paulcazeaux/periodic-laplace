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

#ifndef periodic_laplace_integrate
#define periodic_laplace_integrate

/** \file integrate.hpp Calculation of the boundary integral of grid functions. */

#include "bempp/common/common.hpp"
#include "bempp/common/scalar_traits.hpp"

#include "bempp/assembly/evaluation_options.hpp"
#include "bempp/fiber/quadrature_strategy.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ValueType> class Function;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
class EvaluationOptions;
class GeometryFactory;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
/** \endcond */

/** \relates GridFunction
 *  \brief Calculate the boundary integral of a GridFunction.
 *
 *  This function can be used to estimate the homogenized coefficients from 
 *  a numerical solution to the corrector problem \p trialGridFunction.
 *
 *  The quadrature strategy \p quadStrategy is used to evaluate any
 *  necessary integrals; the \p options object controls the level of
 *  parallelism.
 *
 *  \note If you use the numerical quadrature strategy, you may need to increase
 *  the quadrature order for regular integrals on single elements by at least 2
 *  to ensure that this function produces results with sufficient accuracy.
 */

template <typename BasisFunctionType, typename ResultType>
void Integrate(
        const GridFunction<BasisFunctionType, ResultType>& trialGridFunction,
        const Fiber::Function<ResultType>& testGridFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options,
        typename ScalarTraits<BasisFunctionType>::RealType& integral);

/** \relates GridFunction
 *  \brief Calculate the boundary integral of a GridFunction.
 *
 *  This is an overloaded version, provided for convenience. It is
 *  equivalent to the six-parameter version with \p options set to
 *  <tt>EvaluationOptions()</tt>. */

template <typename BasisFunctionType, typename ResultType>
void Integrate(
        const GridFunction<BasisFunctionType, ResultType>& trialGridFunction,
        const Fiber::Function<ResultType>& testGridFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        typename ScalarTraits<BasisFunctionType>::RealType& integral);

} // namespace Bempp

#endif
