// Copyright (C) 2011-2012 by the BEM++ Authors
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

#include "laplace_3d_periodic_synthetic_boundary_operator_builder.hpp"

#include "bempp/assembly/abstract_boundary_operator.hpp"
#include "bempp/assembly/boundary_operator.hpp"
#include "bempp/assembly/context.hpp"
#include "bempp/assembly/synthetic_nonhypersingular_integral_operator_builder.hpp"

#include "bempp/fiber/explicit_instantiation.hpp"

#include "bempp/space/space.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dPeriodicSyntheticBoundaryOperator(
    BoundaryOperator<BasisFunctionType, ResultType>
    (*constructor)(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& /*context*/,
        const shared_ptr<const Space<BasisFunctionType> >& /*domain*/,
        const shared_ptr<const Space<BasisFunctionType> >& /*range*/,
        const shared_ptr<const Space<BasisFunctionType> >& /*dualToRange*/,
        const std::string& /*label*/,
        int /*symmetry*/),
    const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
    const shared_ptr<const Space<BasisFunctionType> >& domain,
    const shared_ptr<const Space<BasisFunctionType> >& range,
    const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
    std::string label,
    int internalSymmetry,
    int maximumSyntheseSymmetry)
{
    if (!domain || !range || !dualToRange)
        throw std::invalid_argument(
            "makeLaplace3dPeriodicSyntheticBoundaryOperator(): "
            "domain, range and dualToRange must not be null");
    AssemblyOptions internalAssemblyOptions = context->assemblyOptions();
    AcaOptions internalAcaOptions = internalAssemblyOptions.acaOptions();
    internalAcaOptions.mode = AcaOptions::GLOBAL_ASSEMBLY;
    internalAssemblyOptions.switchToAcaMode(internalAcaOptions);
    typedef Context<BasisFunctionType, ResultType> Ctx;
    shared_ptr<Ctx> internalContext(new Ctx(
            context->quadStrategy(), internalAssemblyOptions));
    shared_ptr<const Space<BasisFunctionType> > internalTrialSpace = 
        domain->discontinuousSpace(domain);
    shared_ptr<const Space<BasisFunctionType> > internalTestSpace = 
        dualToRange->discontinuousSpace(dualToRange);

    if (label.empty())
        label = AbstractBoundaryOperator<BasisFunctionType, ResultType>::
            uniqueLabel();

    int syntheseSymmetry = 
        (domain == dualToRange && internalTrialSpace == internalTestSpace) ?
        maximumSyntheseSymmetry : 0;
    // std::cout << "syntheseSymmetry: " << syntheseSymmetry << std::endl;
    BoundaryOperator<BasisFunctionType, ResultType> internalOp =
        constructor(
            internalContext, internalTrialSpace, range /* or whatever */,
            internalTestSpace,
            "(" + label + ")_internal", internalSymmetry);
    return syntheticNonhypersingularIntegralOperator(
            internalOp, domain, range, dualToRange,
            internalTrialSpace, internalTestSpace,
            label, syntheseSymmetry); // TODO: use symmetry too
}

#define INSTANTIATE_FUNCTION(BASIS, RESULT) \
    template BoundaryOperator<BASIS, RESULT> \
    laplace3dPeriodicSyntheticBoundaryOperator( \
        BoundaryOperator<BASIS, RESULT> \
        (*)( \
            const shared_ptr<const Context<BASIS, RESULT> >& /*context*/, \
            const shared_ptr<const Space<BASIS> >& /*domain*/, \
            const shared_ptr<const Space<BASIS> >& /*range*/, \
            const shared_ptr<const Space<BASIS> >& /*dualToRange*/, \
            const std::string& /*label*/, \
            int /*symmetry*/), \
        const shared_ptr<const Context<BASIS, RESULT> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        std::string, int, int)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_FUNCTION);

} // namespace Bempp
