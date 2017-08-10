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

#include "CapsuleBasedDivisionRule.hpp"
#include "RandomNumberGenerator.hpp"
#include "TypeSixSecretionEnumerations.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CapsuleBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // Get separation parameter

    // Make a random direction vector of the required length
    c_vector<double, SPACE_DIM> axis_vector;

    /*
     * Pick a random direction and move the parent cell backwards by 0.5*separation
     * in that direction and return the position of the daughter cell 0.5*separation
     * forwards in that direction.
     */
    switch (SPACE_DIM)
    {
        case 1:
        {
            EXCEPTION("Capsules not coded up in 1D");

        }
        case 2:
        {

        	Node<SPACE_DIM>* pNodeA = rCellPopulation.GetNodeCorrespondingToCell(pParentCell);

        	const double orientation_angle = pNodeA->rGetNodeAttributes()[NA_ANGLE];
        	const double cell_length = pNodeA->rGetNodeAttributes()[NA_LENGTH];


            axis_vector(0) = cell_length*0.25*cos(orientation_angle);
            axis_vector(1) = cell_length*0.25*sin(orientation_angle);
            break;
        }
        case 3:
        {
            EXCEPTION("Capsules not coded up in 3D");

        }
        default:
            // This can't happen
            NEVER_REACHED;
    }

    c_vector<double, SPACE_DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell) - axis_vector;
    c_vector<double, SPACE_DIM> daughter_position = parent_position + 2.0*axis_vector;

    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions(parent_position, daughter_position);

    return positions;
}

// Explicit instantiation
template class CapsuleBasedDivisionRule<1,1>;
template class CapsuleBasedDivisionRule<1,2>;
template class CapsuleBasedDivisionRule<2,2>;
template class CapsuleBasedDivisionRule<1,3>;
template class CapsuleBasedDivisionRule<2,3>;
template class CapsuleBasedDivisionRule<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CapsuleBasedDivisionRule)
