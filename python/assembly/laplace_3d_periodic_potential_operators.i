%{
#define SWIG_FILE_WITH_INIT
#include "assembly/laplace_3d_periodic_single_layer_potential_operator.hpp"
#include "assembly/laplace_3d_periodic_double_layer_potential_operator.hpp"
%}


%feature("compactdefaultargs") laplace3dPeriodicSingleLayerPotentialOperator;
%feature("compactdefaultargs") laplace3dPeriodicDoubleLayerPotentialOperator;

%inline %{

template <typename BasisFunctionType, typename ResultType>
boost::shared_ptr<Bempp::PotentialOperator<BasisFunctionType, ResultType> >
laplace3dPeriodicSingleLayerPotentialOperator()
{
    typedef Bempp::Laplace3dPeriodicSingleLayerPotentialOperator<BasisFunctionType, ResultType> Type;
    return boost::shared_ptr<Type>(new Type);
}

template <typename BasisFunctionType, typename ResultType>
boost::shared_ptr<Bempp::PotentialOperator<BasisFunctionType, ResultType> >
laplace3dPeriodicDoubleLayerPotentialOperator()
{
    typedef Bempp::Laplace3dPeriodicDoubleLayerPotentialOperator<BasisFunctionType, ResultType> Type;
    return boost::shared_ptr<Type>(new Type);
}
%}
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(laplace3dPeriodicSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(laplace3dPeriodicDoubleLayerPotentialOperator);
