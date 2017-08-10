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

#include "ForwardEulerNumericalMethodForCapsules.hpp"
#include "RandomNumberGenerator.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "UblasCustomFunctions.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ForwardEulerNumericalMethodForCapsules<ELEMENT_DIM,SPACE_DIM>::ForwardEulerNumericalMethodForCapsules()
    : AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ForwardEulerNumericalMethodForCapsules<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt)
{
    // Apply forces to each cell, and save a vector of net forces F
    this->ComputeForcesIncludingDamping();

    for (auto node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        double radius = node_iter->rGetNodeAttributes()[NA_RADIUS];
        double length = node_iter->rGetNodeAttributes()[NA_LENGTH];

        double mass = M_PI * radius * radius * (length + 4.0 * radius / 3.0);

        double moment_of_inertia = M_PI * radius * radius * (length * (length * length / 12.0 + radius * radius / 4.0)
                                                             + (4.0 * radius / 3.0) * (0.4 * radius * radius + 0.5 * length * length + 3.0 * length * radius / 8.0));


        node_iter->rGetModifiableLocation() += dt * node_iter->rGetAppliedForce() / mass;
        node_iter->rGetNodeAttributes()[NA_ANGLE] += dt * node_iter->rGetNodeAttributes()[NA_APPLIED_ANGLE] / moment_of_inertia;
    }

//    auto p_rand_gen = RandomNumberGenerator::Instance();

//    for (auto node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
//         node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
//         ++node_iter)
//    {
//        double angle = node_iter->rGetNodeAttributes()[NA_ANGLE];
//
//        // Move a random amount
//        auto random_movement = Create_c_vector(cos(angle), sin(angle));
//        random_movement *= p_rand_gen->NormalRandomDeviate(0.01, 0.005);
//
//        // Rotate a random amount
//        double random_rotation = p_rand_gen->NormalRandomDeviate(0.0, 0.005 * M_PI);
//
//        node_iter->rGetModifiableLocation() += random_movement;
//        node_iter->rGetNodeAttributes()[NA_ANGLE] += random_rotation;
//    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ForwardEulerNumericalMethodForCapsules<ELEMENT_DIM, SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodParameters(rParamsFile);
}

// Explicit instantiation
template class ForwardEulerNumericalMethodForCapsules<1,1>;
template class ForwardEulerNumericalMethodForCapsules<1,2>;
template class ForwardEulerNumericalMethodForCapsules<2,2>;
template class ForwardEulerNumericalMethodForCapsules<1,3>;
template class ForwardEulerNumericalMethodForCapsules<2,3>;
template class ForwardEulerNumericalMethodForCapsules<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ForwardEulerNumericalMethodForCapsules)
