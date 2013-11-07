%{
#define SWIG_FILE_WITH_INIT
#include "assembly/laplace_3d_periodic_single_layer_potential_operator.hpp"
#include "assembly/laplace_3d_periodic_double_layer_potential_operator.hpp"
%}


%feature("compactdefaultargs") Laplace3dPeriodicSingleLayerPotentialOperator;
%feature("compactdefaultargs") Laplace3dPeriodicDoubleLayerPotentialOperator;

#define shared_ptr boost::shared_ptr
%include "assembly/laplace_3d_periodic_single_layer_potential_operator.hpp"
%include "assembly/laplace_3d_periodic_double_layer_potential_operator.hpp"
#undef shared_ptr
namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dPeriodicSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dPeriodicDoubleLayerPotentialOperator);
} // namespace Bempp