/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef _TESTCAPSULEFORCE_HPP_
#define _TESTCAPSULEFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "TypeSixSecretionEnumerations.hpp"

#include "ForwardEulerNumericalMethodForCapsules.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestNumericalMethodForCapsules : public CxxTest::TestSuite
{


public:

    void TestNumericalMethod() throw(Exception)
    {
        ForwardEulerNumericalMethodForCapsules<2, 2> numerical_method;
    }

    void TestCalculateMassOfCapsule() throw(Exception)
    {
        ForwardEulerNumericalMethodForCapsules<2, 2> method;

        // No radius -> no mass
        TS_ASSERT_DELTA(method.CalculateMassOfCapsule(1.23, 0.0), 0.0, 1e-6);

        // No length -> mass of sphere
        TS_ASSERT_DELTA(method.CalculateMassOfCapsule(0.0, 1.23), 7.79478146201, 1e-6);

        // Some odd numbers
        TS_ASSERT_DELTA(method.CalculateMassOfCapsule(2.34, 3.45), 259.506077522, 1e-6);
    }

    void TestCalculateMomentOfInertiaOfCapsule() throw(Exception)
    {
        ForwardEulerNumericalMethodForCapsules<2, 2> method;

        // No radius -> no mass
        TS_ASSERT_DELTA(method.CalculateMomentOfInertiaOfCapsule(1.23, 0.0), 0.0, 1e-6);

        // No length -> moment of inertia of sphere
        TS_ASSERT_DELTA(method.CalculateMomentOfInertiaOfCapsule(0.0, 1.23), 4.71708994955, 1e-6);

        // Calculate MoI for rod length 234, radius 1.23, and rod of length 234 + 2*radius.
        // The real answer is between the two.
        constexpr double r = 1.23;
        constexpr double l = 234.0;
        constexpr double l2 = l + r + r;

        constexpr double lower_bound = (M_PI * r * r * l / 12.0) * (3.0 * r * r + l * l);
        constexpr double upper_bound = (M_PI * r * r * l2 / 12.0) * (3.0 * r * r + l2 * l2);

        TS_ASSERT(method.CalculateMomentOfInertiaOfCapsule(l, r) > lower_bound);
        TS_ASSERT(method.CalculateMomentOfInertiaOfCapsule(l, r) < upper_bound);
    }

};

#endif /*_TESTCAPSULEFORCE_HPP_*/