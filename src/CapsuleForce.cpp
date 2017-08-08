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

#include "CapsuleForce.hpp"

#include "DundeeTypeSixSecretionEnumerations.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point.hpp>

#include "Debug.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CapsuleForce<ELEMENT_DIM,SPACE_DIM>::CapsuleForce()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>()
{
    assert(ELEMENT_DIM == 2u);
    assert(SPACE_DIM == 2u);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CapsuleForce<ELEMENT_DIM,SPACE_DIM>::CalculateDistanceBetweenCapsules(
        Node<SPACE_DIM>& rNodeA,
        Node<SPACE_DIM>& rNodeB)
{
    using geom_point = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;
    using geom_segment = boost::geometry::model::segment<geom_point>;

    auto location_a = rNodeA.rGetLocation();
    auto location_b = rNodeB.rGetLocation();

    const double angle_a = rNodeA.rGetNodeAttributes()[NA_ANGLE];
    const double length_a = rNodeA.rGetNodeAttributes()[NA_LENGTH];

    const double angle_b = rNodeB.rGetNodeAttributes()[NA_ANGLE];
    const double length_b = rNodeB.rGetNodeAttributes()[NA_LENGTH];

    geom_point capsule_a_end_1(location_a[0] + 0.5 * length_a * cos(angle_a),
                               location_a[1] + 0.5 * length_a * sin(angle_a));

    geom_point capsule_a_end_2(location_a[0] - 0.5 * length_a * cos(angle_a),
                               location_a[1] - 0.5 * length_a * sin(angle_a));

    geom_point capsule_b_end_1(location_b[0] + 0.5 * length_b * cos(angle_b),
                               location_b[1] + 0.5 * length_b * sin(angle_b));

    geom_point capsule_b_end_2(location_b[0] - 0.5 * length_b * cos(angle_b),
                               location_b[1] - 0.5 * length_b * sin(angle_b));

    geom_segment capsule_axis_a(capsule_a_end_1, capsule_a_end_2);
    geom_segment capsule_axis_b(capsule_b_end_1, capsule_b_end_2);

    return boost::geometry::distance(capsule_axis_a, capsule_axis_b);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CapsuleForce<ELEMENT_DIM,SPACE_DIM>::CalculateOverlapBetweenCapsules(
        Node<SPACE_DIM>& rNodeA,
        Node<SPACE_DIM>& rNodeB,
        const double distance)
{
    const double radius_a = rNodeA.rGetNodeAttributes()[NA_RADIUS];
    const double radius_b = rNodeB.rGetNodeAttributes()[NA_RADIUS];

    return radius_a + radius_b - distance;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> CapsuleForce<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    c_vector<double, SPACE_DIM> answer;
    return answer;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CapsuleForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class CapsuleForce<1,1>;
template class CapsuleForce<1,2>;
template class CapsuleForce<2,2>;
template class CapsuleForce<1,3>;
template class CapsuleForce<2,3>;
template class CapsuleForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CapsuleForce)
