# Copyright (C) 2013-2014 by Paul Cazeaux
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import bempp.core_periodic_laplace as core_periodic_laplace
from bempp.lib import _constructObjectTemplatedOnBasisAndResult

def _constructOperator(className, context, domain, range, dualToRange, label=None):
    # determine basis function type
    basisFunctionType = domain.basisFunctionType()
    if (basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of all spaces must be the same")

    # determine result type
    resultType = context.resultType()

    if label:
        result = _constructObjectTemplatedOnBasisAndResult(
            core_periodic_laplace, className, basisFunctionType, resultType,
            context, domain, range, dualToRange, label)
    else:
        result = _constructObjectTemplatedOnBasisAndResult(
            core_periodic_laplace, className, basisFunctionType, resultType,
            context, domain, range, dualToRange)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def _constructLaplacePotentialOperator(className, context):
    basisFunctionType = context.basisFunctionType()
    resultType = context.resultType()
    result = _constructObjectTemplatedOnBasisAndResult(
        core_periodic_laplace, className, basisFunctionType, resultType)
    result._context = context
    return result

# determine the type used to represent the values of the basis functions into
# which functions acted upon by the operator will be expanded, and the type used
# to represent the values of the functions produced by this operator.

def createLaplace3dPeriodicSingleLayerBoundaryOperator(
                                               context, domain, range, dualToRange, label=None):
    """
        Create and return a single-layer-potential boundary operator for the
        Laplace equation in a periodic box in 3D.
        
        *Parameters:*
        - context (Context)
        A Context object to control the assembly of the weak form of the
        newly constructed operator.
        - domain (Space)
        Function space to be taken as the domain of the operator.
        - range (Space)
        Function space to be taken as the range of the operator.
        - dualToRange (Space)
        Function space to be taken as the dual to the range of the operator.
        - label (string)
        Textual label of the operator. If set to None (default), a unique
        label will be generated automatically.
        
        *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
        object, with BasisFunctionType and ResultType determined automatically from
        the context argument and equal to either float32, float64, complex64 or
        complex128.
        """
    return _constructOperator(
                              "laplace3dPeriodicSingleLayerBoundaryOperator",
                              context, domain, range, dualToRange, label)

def createLaplace3dPeriodicDoubleLayerBoundaryOperator(
                                               context, domain, range, dualToRange, label=None):
    """
        Create and return a double-layer-potential boundary operator for the
        Laplace equation in a periodic box in 3D.
        
        *Parameters:*
        - context (Context)
        A Context object to control the assembly of the weak form of the
        newly constructed operator.
        - domain (Space)
        Function space to be taken as the domain of the operator.
        - range (Space)
        Function space to be taken as the range of the operator.
        - dualToRange (Space)
        Function space to be taken as the dual to the range of the operator.
        - label (string)
        Textual label of the operator. If set to None (default), a unique
        label will be generated automatically.
        
        *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
        object, with BasisFunctionType and ResultType determined automatically from
        the context argument and equal to either float32, float64, complex64 or
        complex128.
        """
    return _constructOperator(
                              "laplace3dPeriodicDoubleLayerBoundaryOperator",
                              context, domain, range, dualToRange, label)

def createLaplace3dPeriodicAdjointDoubleLayerBoundaryOperator(
                                                      context, domain, range, dualToRange, label=None):
    """
        Create and return an adjoint double-layer-potential boundary operator for
        the Laplace equation in a periodic box in 3D.
        
        *Parameters:*
        - context (Context)
        A Context object to control the assembly of the weak form of the
        newly constructed operator.
        - domain (Space)
        Function space to be taken as the domain of the operator.
        - range (Space)
        Function space to be taken as the range of the operator.
        - dualToRange (Space)
        Function space to be taken as the dual to the range of the operator.
        - label (string)
        Textual label of the operator. If set to None (default), a unique
        label will be generated automatically.
        
        *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
        object, with BasisFunctionType and ResultType determined automatically from
        the context argument and equal to either float32, float64, complex64 or
        complex128.
        """
    return _constructOperator(
                              "laplace3dPeriodicAdjointDoubleLayerBoundaryOperator",
                              context, domain, range, dualToRange, label)

def createLaplace3dPeriodicSingleLayerPotentialOperator(context):
    """
        Create and return a single-layer potential operator for the Laplace
        equation in a periodic box in 3D.
        
        *Parameters:*
        - context (Context)
        A Context object used to control the evaluation of integrals
        occurring in the definition of the potential operator.
        
        *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
        object, with BasisFunctionType and ResultType determined automatically from
        the context argument and equal to either float32, float64, complex64 or
        complex128.
        
        Note about BEM++ terminology: a *potential operator* acts on functions
        defined on a surface S and produces functions defined at any point of the
        space surrounding S, but not necessarily on S itself. In contrast, a
        *boundary operator* acts on on functions defined on a surface S and produces
        functions defined on the same surface S.
        """
    return _constructLaplacePotentialOperator(
                                              "laplace3dPeriodicSingleLayerPotentialOperator", context)

def createLaplace3dPeriodicDoubleLayerPotentialOperator(context):
    """
        Create and return a double-layer potential operator for the Laplace
        equation in a periodic box in 3D.
        
        *Parameters:*
        - context (Context)
        A Context object used to control the evaluation of integrals
        occurring in the definition of the potential operator.
        
        *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
        object, with BasisFunctionType and ResultType determined automatically from
        the context argument and equal to either float32, float64, complex64 or
        complex128.
        
        Note about BEM++ terminology: a *potential operator* acts on functions
        defined on a surface S and produces functions defined at any point of the
        space surrounding S, but not necessarily on S itself. In contrast, a
        *boundary operator* acts on on functions defined on a surface S and produces
        functions defined on the same surface S.
        """
    return _constructLaplacePotentialOperator(
                                              "laplace3dPeriodicDoubleLayerPotentialOperator", context)

    def Integrate(
        gridFunction, testFunction, quadStrategy, evaluationOptions,
        surfaceNormalDependent=False):
    """
    Calculate the integral of a product between a grid function and a test function
    defined as a Python callable.

    *Parameters:*
       - gridFunction (GridFunction)
            A grid function.
       - testFunction (a Python callable object)
            A Python callable object whose values on 'gridFunction.grid()' will
            be multiplied by those of the grid function. If
            'surfaceNormalDependent' is set to False (default), 'testFunction'
            will be passed a single argument containing a 1D array of the
            coordinates of a point lying on the grid 'gridFunction.grid()'. If
            'surfaceNormalDependent' is set to True, the function will be passed
            one more argument, a 1D array containing the components of the unit
            vector normal to 'gridFunction.grid()' at the point given in the
            first argument. In both cases 'testFunction' should return its
            value at the given point, in the form of a scalar or a 1D array with
            dimension equal to 'space.codomainDimension()'.
       - quadStrategy (QuadratureStrategy)
            Quadrature stategy to be used to evaluate integrals involved in the
            boundary integral.
       - evaluationOptions (EvaluationOptions)
            Additional options controlling function evaluation.
       - surfaceNormalDependent (bool)
            Indicates whether 'testFunction' depends on the unit vector
            normal to the grid or not.

    *Returns* the value of the boundary integralof the product gridFunction * testFunction.

    See the documentation of createGridFunction() for example definitions of
    Python functions that can be passed to the 'testFunction' argument.
    """
    basisFunctionType = checkType(gridFunction.basisFunctionType())
    resultType = checkType(gridFunction.resultType())
    if (basisFunctionType != quadStrategy.basisFunctionType() or
            resultType != quadStrategy.resultType()):
        raise TypeError("BasisFunctionType and ResultType of gridFunction "
                        "and testFunction must be the same")
    if surfaceNormalDependent:
        dependent = "Dependent"
    else:
        dependent = "Independent"
    functor = _constructObjectTemplatedOnValue(
        core, "PythonSurfaceNormal%sFunctor" % dependent,
        resultType, testFunction,
        gridFunction.space().grid().dimWorld(), # argument dimension
        gridFunction.space().codomainDimension() # result dimension
        )
    return _constructObjectTemplatedOnBasisAndResult(
        core, "IntegralFromPythonSurfaceNormal%sFunctor" % dependent,
        basisFunctionType, resultType,
        gridFunction, functor, quadStrategy, evaluationOptions)
