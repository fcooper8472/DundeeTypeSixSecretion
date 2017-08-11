/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "SquareBoundaryCondition.hpp"

#include "CellLabel.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "Debug.hpp"

SquareBoundaryCondition::SquareBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation)
        : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation)
{
    if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation) == NULL)
    {
        EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
    }
}

void SquareBoundaryCondition::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    // Iterate over the cell population
    for (auto cell_iter = this->mpCellPopulation->Begin(); cell_iter != this->mpCellPopulation->End(); ++cell_iter)
    {
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        Node<2>* p_node = this->mpCellPopulation->GetNode(node_index);
        c_vector<double,2>& cell_location = p_node->rGetModifiableLocation();

        if (cell_location[0] < 0.0)
        {
            cell_location[0] = 0.0;
        }
        if ((cell_location[0] > 10.0))
        {
            cell_location[0] = 10.0;
        }
        if ((cell_location[1] < 0.0))
        {
            cell_location[1] = 0.0;
        }
        if ((cell_location[1] > 10.0))
        {
            cell_location[1] = 10.0;
        }
    }
}

bool SquareBoundaryCondition::VerifyBoundaryCondition()
{
    return true;
}

void SquareBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SquareBoundaryCondition)
