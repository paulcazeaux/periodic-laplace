%{
#define SWIG_FILE_WITH_INIT
#include "assembly/integrate.hpp"
%}

%inline %{

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
typename ScalarTraits<BasisFunctionType>::RealType
IntegrateFromPythonSurfaceNormalIndependentFunctor(
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
typename ScalarTraits<BasisFunctionType>::RealType
IntegrateFromPythonSurfaceNormalDependentFunctor(
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