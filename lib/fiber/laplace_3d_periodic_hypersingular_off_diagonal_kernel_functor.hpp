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

#ifndef fiber_laplace_3d_periodic_hypersingular_off_diagonal_kernel_functor_hpp
#define fiber_laplace_3d_periodic_hypersingular_off_diagonal_kernel_functor_hpp

#include "bempp/common/common.hpp"

#include "bempp/fiber/geometrical_data.hpp"
#include "bempp/fiber/scalar_traits.hpp"
#include "boost/math/special_functions/spherical_harmonic.hpp"

namespace Fiber
{

/** \ingroup laplace_3d_periodic
 *  \ingroup functors
 *  \brief Kernel functor for the hypersingular operator associated with the
 *  Laplace equation in a periodic box of size 1 in 3D applicable for test and trial points
 *  lying on nonadjacent elements.
 *  This functor computes the regular part of the kernel corresponding to the 26 Green's function
 *  in the neighboring cells and a correction part approximated using tabulated spherical harmonics.
 *  
 *
 *  \tparam ValueType Type used to represent the values of the kernel. It can
 *  be one of: \c float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>.
 *
 *  \see laplace_3d_periodic
 */

template <typename ValueType_>
class Laplace3dPeriodicHypersingularOffDiagonalKernelFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int kernelCount() const { return 1; }
    int kernelRowCount(int /* kernelIndex */) const { return 1; }
    int kernelColCount(int /* kernelIndex */) const { return 1; }

    void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const {
        testGeomDeps |= GLOBALS | NORMALS;
        trialGeomDeps |= GLOBALS | NORMALS;
    }

    template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
    void evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& testGeomData,
            const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            CollectionOf2dSlicesOfNdArrays<ValueType>& result) const {
        const int coordCount = 3;
        assert(testGeomData.dimWorld() == coordCount);
        assert(result.size() == 1);
        
        /*  Define coefficients for the regular part of the periodic Green's kernel for the 3D Laplace equation in a periodic box of size 1.
         *      We use here an expansion with harmonics up to degree 20, which yields an approximation correct up to precision 1e-10.
         *      The singular part contains the 27 Green's function for the neighbouring cells and a corrective term for a uniform neutralizing
         *      background charge.
         */
        
        const CoordinateType Harmonic_Coefficients[34] = {
            0.01101709645861367,
            0.01316794887913183,
            5.954412827470426e-05,
            -0.0002227937273979395,
            0.0002062746586075765,
            0.0001551394663489125,
            0.0002363741367920927,
            5.288864558642977e-06,
            -1.065878569810829e-05,
            -1.268649253304447e-05,
            2.855411046592645e-05,
            2.503036078346336e-05,
            1.043475591993333e-05,
            3.395059459839415e-05,
            1.255601659344675e-06,
            -1.846641841052362e-06,
            -1.981606873099481e-06,
            -2.406130170392715e-06,
            1.25747200483749e-06,
            6.475400668217489e-07,
            7.957080434716371e-07,
            9.28658402338717e-07,
            1.321949011147113e-06,
            9.74822898071474e-08,
            -1.330176890868469e-07,
            -4.957180138283109e-08,
            -2.127820763047173e-07,
            -1.370001191920438e-07,
            2.762838013591894e-08,
            -5.682089767011293e-09,
            1.844902148801811e-08,
            2.90202283948248e-08,
            3.429155607637934e-08,
            1.998351880226059e-08

        };
        
        ValueType ri2[3*coordCount];
        ValueType res = 0;
        CoordinateType numeratorCoeff[3*coordCount];
        CoordinateType diff[3];
        
        const CoordinateType common = static_cast<CoordinateType>(-1. / (4. * M_PI));
        const CoordinateType ONE = static_cast<CoordinateType>(1.).
        CoordinateType nTest_diff[3*coordCount];
        CoordinateType nTrial_diff[3*coordCount];
        CoordinateType nTest_nTrial = 0.;
        
        for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
            diff[coordIndex] = boost::math::tools::fmod_workaround(testGeomData.global(coordIndex)
                                                                   -  trialGeomData.global(coordIndex)
                                                                   + static_cast<CoordinateType>(0.5), ONE);
            if (diff[coordIndex] > 0)
                diff[coordIndex] -=  static_cast<CoordinateType>(0.5);
            else
                diff[coordIndex] +=  static_cast<CoordinateType>(0.5);    // Periodization
            
            
            ri2[  3*coordIndex] = (diff[coordIndex]-ONE)*(diff[coordIndex]-ONE);
            ri2[1+3*coordIndex] = diff[coordIndex]*diff[coordIndex];
            ri2[2+3*coordIndex] = (diff[coordIndex]+ONE)*(diff[coordIndex]+ONE);
            
            nTest_diff[  3*coordIndex] = (diff[coordIndex]-ONE) * testGeomData.normal(coordIndex);
            nTest_diff[1+3*coordIndex] = diff[coordIndex] * testGeomData.normal(coordIndex);
            nTest_diff[2+3*coordIndex] = (diff[coordIndex]+ONE) * testGeomData.normal(coordIndex);
            nTrial_diff[  3*coordIndex] = (diff[coordIndex]-ONE) * trialGeomData.normal(coordIndex);
            nTrial_diff[1+3*coordIndex] = diff[coordIndex] * trialGeomData.normal(coordIndex);
            nTrial_diff[2+3*coordIndex] = (diff[coordIndex]+ONE) * trialGeomData.normal(coordIndex);
            
            nTest_nTrial += testGeomData.normal(coordIndex) *
            trialGeomData.normal(coordIndex);
        }
        // Loop over the 27 adjacent cells (!)
        CoordinateType distance2, distance5, numerator;
        const CoordinateType THREE = static_cast<CoordinateType>(3.);
        
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 3; j < 6; ++j)
            {
                for (int k = 6; k < 9; ++k)
                {
                    distance2 = ri2[i] + ri2[j] + ri2[k];
                    distance5 = distance2 * distance2 * sqrt(distance2);
                    numerator = nTest_nTrial * distance2 - THREE * (nTest_diff[i] + nTest_diff[j] + nTest_diff[k]) *
                    (nTrial_diff[i] + nTrial_diff[j] + nTrial_diff[k]));
                }
            }
        }
        
        // corrective term for a neutral background charge
        
        res += static_cast<CoordinateType>(1./3.) * nTest_nTrial;
        
        
        // Approximation of the correction using spherical harmonics
        int ind = 0;
        // Transform to spherical coordinates
        CoordinateType rfactor = sqrt(r2);
        CoordinateType theta = acos(diff[3]/sqrt(r2));
        CoordinateType phi = atan2(diff[2], diff[1]);
        
        for (int degree = 4; degree < 21; degree += 2)
        {
            rfactor = rfactor * r2;
            for (int order = 0; order < degree + 1; order += 4)
            {
                // Compute the spherical harmonics we need in the gradient calculation
                std::complex<double> RPlus = static_cast<CoordinateType>(sqrt((2*l+1)*(l-m)*(l-m-1)/(l-1))) *
                boost::math::spherical_harmonic(degree-1, order+1, theta, phi);
                std::complex<double> RMinus = static_cast<CoordinateType>(sqrt((2*l+1)*(l+m)*(l+m-1)/(l-1))) *
                boost::math::spherical_harmonic(degree-1, order-1, theta, phi);
                CoordinateType Rz =  static_cast<CoordinateType>(sqrt((2*l+1)*(l+m)*(l-m)/(2*l-1))) *
                boost::math::spherical_harmonic_r(degree-1, order, theta, phi);
                
                // x component
                res += Harmonic_Coefficients[ind] * (RPlus.real() - RMinus.real()) * testGeomData.normal(0) * rfactor;
                // y component
                res -= Harmonic_Coefficients[ind] * (RPlus.imag() + RMinus.imag()) * testGeomData.normal(1) * rfactor;
                // z component
                res += Harmonic_Coefficients[ind] * Rz * testGeomData.normal(2) * rfactor;
                ++ind;
            }
        }

        
        result[0](0, 0) = res;

//        CoordinateType distanceSq = 0.;
//        CoordinateType nTest_diff = 0;
//        CoordinateType nTrial_diff = 0;
//        CoordinateType nTest_nTrial = 0.;
//        for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
//            CoordinateType diff = trialGeomData.global(coordIndex) -
//                    testGeomData.global(coordIndex);
//            distanceSq += diff * diff;
//            nTest_diff += diff * testGeomData.normal(coordIndex);
//            nTrial_diff += diff * trialGeomData.normal(coordIndex);
//            nTest_nTrial += testGeomData.normal(coordIndex) *
//            trialGeomData.normal(coordIndex);
//        }
//        CoordinateType distance = sqrt(distanceSq);
//        CoordinateType commonFactor =
//                static_cast<CoordinateType>(-1. / (4. * M_PI)) /
//                (distance * distanceSq * distanceSq);
//        result[0](0, 0) = commonFactor * (
//                    nTest_nTrial * distanceSq -
//                    static_cast<CoordinateType>(3.) * nTest_diff * nTrial_diff);
    }
};

} // namespace Fiber

#endif
