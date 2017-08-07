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

#include "DundeeTypeSixSecretionEnumerations.hpp"

#include "OutputFileHandler.hpp"
#include "Node.hpp"
#include "NodeAttributes.hpp"
#include "PetscTools.hpp"

#include "CapsuleForce.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCapsuleForce : public CxxTest::TestSuite
{


public:

    void TestDistanceBetweenTwoCapsules() throw(Exception)
    {
        CapsuleForce<2, 2> force;

        // Two horizontal rods a distance 2 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{1.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_ANGLE] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<2> node_b(0u, std::vector<double>{1.0, 2.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_ANGLE] = 0.0;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.25;

            TS_ASSERT_DELTA(force.CalculateDistanceBetweenCapsules(node_a, node_b), 1.25, 1e-6);
        }

        // One horizontal, one vertical, distance 2 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{1.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_ANGLE] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.75;

            Node<2> node_b(0u, std::vector<double>{1.0, 3.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_ANGLE] = 0.5 * M_PI;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.05;

            TS_ASSERT_DELTA(force.CalculateDistanceBetweenCapsules(node_a, node_b), 1.2, 1e-6);
        }

        // Intersecting cross shape, distance 0 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{0.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_ANGLE] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<2> node_b(0u, std::vector<double>{0.0, 0.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_ANGLE] = 0.5 * M_PI;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.5;

            TS_ASSERT_DELTA(force.CalculateDistanceBetweenCapsules(node_a, node_b), -1.0, 1e-6);
        }

        // Intersecting cross shape, distance 0 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{4.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_ANGLE] = 0.0;
            attributes_a[NA_LENGTH] = 4.0;
            attributes_a[NA_RADIUS] = 0.12;

            Node<2> node_b(0u, std::vector<double>{4.0, 3.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_ANGLE] = atan(0.75);
            attributes_b[NA_LENGTH] = 4.0;
            attributes_b[NA_RADIUS] = 0.13;

            TS_ASSERT_DELTA(force.CalculateDistanceBetweenCapsules(node_a, node_b), 1.55, 1e-6);
        }
    }
};

#endif /*_TESTCAPSULEFORCE_HPP_*/