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
from bempp.lib import _constructOperator, _constructLaplacePotentialOperator

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
