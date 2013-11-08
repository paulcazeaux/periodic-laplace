%{
#define SWIG_FILE_WITH_INIT
#include "assembly/laplace_3d_periodic_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_periodic_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_periodic_adjoint_double_layer_boundary_operator.hpp"
%}
%include "numpy.i"
%include "bempp.swg"

%feature("compactdefaultargs") laplace3dPeriodicSingleLayerBoundaryOperator;
%feature("compactdefaultargs") laplace3dPeriodicDoubleLayerBoundaryOperator;
%feature("compactdefaultargs") laplace3dPeriodicAdjointDoubleLayerBoundaryOperator;


#define shared_ptr boost::shared_ptr
%include "assembly/laplace_3d_periodic_single_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_periodic_double_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_periodic_adjoint_double_layer_boundary_operator.hpp"
#undef shared_ptr
namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(laplace3dPeriodicSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(laplace3dPeriodicDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(laplace3dPeriodicAdjointDoubleLayerBoundaryOperator);
} // namespace Bempp