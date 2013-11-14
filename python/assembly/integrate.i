%{
#include "assembly/integrate.hpp"
#include "assembly/surface_normal_dependent_function.hpp"
#include "assembly/surface_normal_independent_function.hpp"
%}

%inline %{

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
ResultType IntegrateFromPythonSurfaceNormalIndependentFunctor(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        const PythonSurfaceNormalIndependentFunctor<ResultType>& functor,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options)
{
    return Integrate(gridFunction,
                              surfaceNormalIndependentFunction(functor),
                              quadStrategy,
                              options);
}
template <typename BasisFunctionType, typename ResultType>
ResultType IntegrateFromPythonSurfaceNormalDependentFunctor(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        const PythonSurfaceNormalDependentFunctor<ResultType>& functor,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options)
{
    return Integrate(gridFunction,
                              surfaceNormalDependentFunction(functor),
                              quadStrategy,
                              options);
}

} // namespace Bempp

%}

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(IntegrateFromPythonSurfaceNormalIndependentFunctor);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(IntegrateFromPythonSurfaceNormalDependentFunctor);

} // namespace Bempp