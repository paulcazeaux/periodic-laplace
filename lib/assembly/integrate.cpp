// Copyright (C) 2013-2014 by the Paul Cazeaux
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

#include "integrate.hpp"

#include "bempp/assembly/evaluation_options.hpp"
#include "bempp/assembly/grid_function.hpp"
#include "bempp/assembly/local_assembler_construction_helper.hpp"

#include "bempp/common/armadillo_fwd.hpp"
#include "bempp/common/complex_aux.hpp"
#include "bempp/common/shared_ptr.hpp"

#include "bempp/fiber/collection_of_2d_arrays.hpp"
#include "bempp/fiber/collection_of_4d_arrays.hpp"
#include "bempp/fiber/default_collection_of_basis_transformations.hpp"
#include "bempp/fiber/default_collection_of_kernels.hpp"
#include "bempp/fiber/default_kernel_trial_integral.hpp"
#include "bempp/fiber/evaluator_for_integral_operators.hpp"
#include "bempp/fiber/explicit_instantiation.hpp"
#include "bempp/fiber/function.hpp"
#include "bempp/fiber/geometrical_data.hpp"

namespace Bempp
{

namespace
{

template <typename ValueType_>
class KernelFunctorFromTestFunction
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    explicit KernelFunctorFromTestFunction(
        const Fiber::Function<ValueType>& function) : m_function(function)
    {}

    int kernelCount() const {
        return 1;
    }

    int kernelRowCount(int kernelIndex) const {
        return m_function.codomainDimension();
    }

    int kernelColCount(int kernelIndex) const {
        return 1;
    }

    void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const {
        m_function.addGeometricalDependencies(trialGeomDeps);
    }

    template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
    void evaluate(
            const Fiber::ConstGeometricalDataSlice<CoordinateType>& testGeomData,
            const Fiber::ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            CollectionOf2dSlicesOfNdArrays<ValueType>& result) const {
        Fiber::GeometricalData<CoordinateType> geomData =
                trialGeomData.asGeometricalData();
        arma::Mat<ValueType> resultMatrix;
        m_function.evaluate(geomData, resultMatrix);
        assert(resultMatrix.n_rows == result[0].extent(0));
        assert(resultMatrix.n_cols == result[0].extent(1));
        for (size_t j = 0; j < resultMatrix.n_cols; ++j)
            for (size_t i = 0; i < resultMatrix.n_rows; ++i)
                result[0](i, j) = resultMatrix(i, j);
    }

private:
    const Fiber::Function<ValueType>& m_function;
};

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class IntegrandFunctor
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    void addGeometricalDependencies(size_t& trialGeomDeps) const {
        // do nothing
    }

    int resultDimension() const {
        return 1;
    }

    template <typename CollectionOf1dSlicesOfConstNdArrays>
    void evaluate(
            const Fiber::ConstGeometricalDataSlice<CoordinateType>& /* trialGeomData */,
            const Fiber::CollectionOf2dSlicesOfConst4dArrays<KernelType>& kernelValues,
            const CollectionOf1dSlicesOfConstNdArrays& trialValues,
            std::vector<ResultType>& result) const {
        // Assert that there is only a single kernel with a single column
        assert(kernelValues.size() == 1);
        assert(kernelValues[0].extent(1) == 1);
        const size_t componentCount = kernelValues[0].extent(0);

        // Assert that there is only one trial transformation
        // and that its dimensions are correct
        assert(trialValues.size() == 1);
#ifndef NDEBUG
        assert(trialValues[0].extent(0) == componentCount);
#endif

        // Assert that the result has only one element
        assert(result.size() == 1);

        // (t - k)* . (t - k)
        result[0] = 0.;
        for (size_t i = 0; i < componentCount; ++i) {
            result[0] += trialValues[0](i) * kernelValues[0](i, 0);
        }
    }
};

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<Fiber::EvaluatorForIntegralOperators<ResultType> >
makeEvaluator(
        const GridFunction<BasisFunctionType, ResultType>& trialGridFunction,
        const Fiber::Function<ResultType>& testFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options)
{
    // Collect the standard set of data necessary for construction of
    // evaluators and assemblers
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Shapeset<BasisFunctionType>*> ShapesetPtrVector;
    typedef std::vector<std::vector<ResultType> > CoefficientsVector;
    typedef LocalAssemblerConstructionHelper Helper;

    shared_ptr<RawGridGeometry> rawGeometry;
    shared_ptr<GeometryFactory> geometryFactory;
    shared_ptr<Fiber::OpenClHandler> openClHandler;
    shared_ptr<ShapesetPtrVector> shapesets;

    const Space<BasisFunctionType>& space = *trialGridFunction.space();
    Helper::collectGridData(space,
                            rawGeometry, geometryFactory);
    Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
                              rawGeometry, openClHandler);
    Helper::collectShapesets(space, shapesets);

    // In addition, get coefficients of argument's expansion in each element
    const GridView& view = space.gridView();
    const int elementCount = view.entityCount(0);

    shared_ptr<CoefficientsVector> localCoefficients =
            boost::make_shared<CoefficientsVector>(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();
    for (int i = 0; i < elementCount; ++i) {
        const Entity<0>& element = it->entity();
        trialGridFunction.getLocalCoefficients(element, (*localCoefficients)[i]);
        it->next();
    }

    // Construct kernels collection
    typedef KernelFunctorFromTestFunction<ResultType> KernelFunctor;
    typedef Fiber::DefaultCollectionOfKernels<KernelFunctor> Kernels;
    shared_ptr<Fiber::CollectionOfKernels<ResultType> > kernels =
        boost::make_shared<Kernels>(KernelFunctor(testFunction));

    // Construct trial function transformations collection
    const Fiber::CollectionOfShapesetTransformations<CoordinateType>&
        trialTransformations = space.basisFunctionValue();

    // Construct integrand
    typedef IntegrandFunctor<BasisFunctionType, ResultType, ResultType>
            IntegralFunctor;
    shared_ptr<Fiber::KernelTrialIntegral<
    BasisFunctionType, ResultType, ResultType> > integrand =
        boost::make_shared<Fiber::DefaultKernelTrialIntegral<IntegralFunctor> >(
            IntegralFunctor());

    // Now create the evaluator
    return quadStrategy.makeEvaluatorForIntegralOperators(
                geometryFactory, rawGeometry,
                shapesets,
                kernels,
                make_shared_from_ref(trialTransformations), // lives as long as space does
                integrand,
                localCoefficients,
                openClHandler,
                options.parallelizationOptions());
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
ResultType Integrate(
        const GridFunction<BasisFunctionType, ResultType>& trialGridFunction,
        const Fiber::Function<ResultType>& testFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options)
{
    // First, construct the evaluator.

    typedef Fiber::EvaluatorForIntegralOperators<ResultType> Evaluator;
    std::auto_ptr<Evaluator> evaluator =
        makeEvaluator(
            trialGridFunction, testFunction, quadStrategy, options);

    typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
    arma::Mat<CoordinateType> evaluationPoints(1, 1);
    evaluationPoints(0, 0) = 0.;
    arma::Mat<ResultType> resultMatrix;

    evaluator->evaluate(Evaluator::FAR_FIELD, evaluationPoints, resultMatrix);

    ResultType result = resultMatrix(0, 0);
    return result;
}

template <typename BasisFunctionType, typename ResultType>
void Integrator(
        const GridFunction<BasisFunctionType, ResultType>& trialGridFunction,
        const Fiber::Function<ResultType>& testFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options,
        ResultType& integral)
{
    integral = Integrate(trialGridFunction, testFunction,
                                                quadStrategy, options);
}

template <typename BasisFunctionType, typename ResultType>
void Integrator(
        const GridFunction<BasisFunctionType, ResultType>& trialGridFunction,
        const Fiber::Function<ResultType>& testFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        ResultType& integral)
{
    Integrator(trialGridFunction, testFunction, quadStrategy,
                    EvaluationOptions(), integral);
}

#define INSTANTIATE_FUNCTION(BASIS, RESULT) \
    template \
    RESULT Integrate( \
            const GridFunction<BASIS, RESULT>& trialGridFunction, \
            const Fiber::Function<RESULT>& testFunction, \
            const Fiber::QuadratureStrategy< \
                BASIS, RESULT, GeometryFactory>& quadStrategy, \
            const EvaluationOptions& options); \
    template \
    void Integrator( \
            const GridFunction<BASIS, RESULT>& trialGridFunction, \
            const Fiber::Function<RESULT>& testFunction, \
            const Fiber::QuadratureStrategy< \
                BASIS, RESULT, GeometryFactory>& quadStrategy, \
            const EvaluationOptions& options, \
            RESULT& integral); \
    template \
    void Integrator( \
            const GridFunction<BASIS, RESULT>& trialGridFunction, \
            const Fiber::Function<RESULT>& testFunction, \
            const Fiber::QuadratureStrategy< \
                BASIS, RESULT, GeometryFactory>& quadStrategy, \
            RESULT& integral)

FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_FUNCTION);

} // namespace Bempp
