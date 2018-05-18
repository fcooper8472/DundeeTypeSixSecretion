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
#include "NodeBasedCellPopulationWithCapsules.hpp"
#include "ReplicatableVector.hpp"
#include "OdeLinearSystemSolver.hpp"
#include "UniformCellCycleModel.hpp"
#include "TypeSixSecretionEnumerations.hpp"

#include "Debug.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"

#include "AbstractCentreBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "TypeSixMachineProperty.hpp"

template<unsigned DIM>
NodeBasedCellPopulationWithCapsules<DIM>::NodeBasedCellPopulationWithCapsules(NodesOnlyMesh<DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh)
    : NodeBasedCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh)
{

}

template<unsigned DIM>
NodeBasedCellPopulationWithCapsules<DIM>::NodeBasedCellPopulationWithCapsules(NodesOnlyMesh<DIM>& rMesh)
    : NodeBasedCellPopulation<DIM>(rMesh)
{
    // No Validate() because the cells are not associated with the cell population yet in archiving
}

//template<unsigned DIM>
//void NodeBasedCellPopulationWithCapsules<DIM>::UpdateNodeLocations(double dt)
//{
//	MARK;
//
//	// Update cell lengths
//    // Iterate over cell population
//
//    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
//       cell_iter != this->End();
//       ++cell_iter)
//    {
//
//
//          double cell_age  = cell_iter->GetAge();//+SimulationTime::Instance()->GetTimeStep();
//          double initial_length = 2.0;
//
//          if (bool(dynamic_cast<UniformCellCycleModel*>(cell_iter->GetCellCycleModel())))
//          {
//              double cell_cycle_time = static_cast<UniformCellCycleModel*>(cell_iter->GetCellCycleModel())->GetCellCycleDuration();
//
//              //if (cell_age >= cell_cycle_time)
//              //{
//            //	  cell_age=0.0;
//            //	  std::cout<< "Ahh " << "\n";
//             // }
//
//
//              //assert(bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&mrCellPopulation)));
//              //Node<DIM>* pNodeA = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)->GetNodeCorrespondingToCell(*cell_iter);
//
//              Node<DIM>* pNodeA = this->GetNodeCorrespondingToCell(*cell_iter);
//              double division_length = 2*initial_length + 2*pNodeA->rGetNodeAttributes()[NA_RADIUS];
//              double new_length = initial_length + (division_length - initial_length)*cell_age/cell_cycle_time;
//              //double new_length = initial_length*(1.0+cell_age/cell_cycle_time);
//
//
//              pNodeA->rGetNodeAttributes()[NA_LENGTH] = new_length;
//
//
//          }
//    }
//
//
//    NodeBasedCellPopulation<DIM>::UpdateNodeLocations(dt);
//}

template<unsigned DIM>
CellPtr NodeBasedCellPopulationWithCapsules<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{

	auto pNewCellTemp=NodeBasedCellPopulation<DIM>::AddCell(pNewCell, pParentCell);


	// Get new node
	Node<DIM>* p_new_node = this->GetNodeCorrespondingToCell(pNewCellTemp);// new Node<DIM>(this->GetNumNodes(), daughter_position, false); // never on boundary

	p_new_node->AddNodeAttribute(0.0);
	p_new_node->rGetNodeAttributes().resize(NA_VEC_LENGTH);

	double angle = (this->GetNodeCorrespondingToCell(pParentCell))->rGetNodeAttributes()[NA_THETA];
	angle = angle + 0.001*(RandomNumberGenerator::Instance()->ranf()-0.5)*2*M_PI;

	double length = (this->GetNodeCorrespondingToCell(pParentCell))->rGetNodeAttributes()[NA_LENGTH];
	double radius = (this->GetNodeCorrespondingToCell(pParentCell))->rGetNodeAttributes()[NA_RADIUS];

	p_new_node->rGetNodeAttributes()[NA_THETA] =  angle;
	if (DIM==3)
	{
		double phi = (this->GetNodeCorrespondingToCell(pParentCell))->rGetNodeAttributes()[NA_PHI];
		p_new_node->rGetNodeAttributes()[NA_PHI] =  phi;
	}

	p_new_node->rGetNodeAttributes()[NA_LENGTH] = length;
	p_new_node->rGetNodeAttributes()[NA_RADIUS] = radius;

    // Get this cell's type six machine property data
    CellPropertyCollection collection = pParentCell->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
    if (collection.GetSize() != 1)
	{
		EXCEPTION("TypeSixMachineCellKiller cannot be used unless each cell has a TypeSixMachineProperty");
	}
	boost::shared_ptr<TypeSixMachineProperty> p_parent_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
	std::vector<std::pair<unsigned, std::vector<double>> >& r_parent_data = p_parent_property->rGetMachineData();



	MAKE_PTR(TypeSixMachineProperty, p_new_daughter_property);
	MAKE_PTR(TypeSixMachineProperty, p_new_parent_property);


	unsigned parent_node_index = this->GetLocationIndexUsingCell(pParentCell);
	Node<DIM>* p_parent_node = this->GetNode(parent_node_index);
	double L = p_parent_node->rGetNodeAttributes()[NA_LENGTH];

	// Iterate over machines in this cell and distribute to mother or daughter
	for (auto& r_pair : r_parent_data)
	{
	        // retrieve machine coordinates in frame of old cell
			std::vector<double> local_machine_coords = r_pair.second;


			double vertical_coordinate=local_machine_coords[0];

			if (vertical_coordinate < 0.0 ) // machine inherited by daughter cell
			{
				double new_vertical_coord = vertical_coordinate+L/4.0;
                double new_phi =  local_machine_coords[1];
                std::vector<double> new_machine_coordinates;
                new_machine_coordinates.push_back(new_vertical_coord);
                new_machine_coordinates.push_back(new_phi);
                p_new_daughter_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(r_pair.first, new_machine_coordinates));
			}
			else if (vertical_coordinate >=0.0)// machine inherited by mother cell
			{

				double new_vertical_coord = vertical_coordinate-L/4.0;
				double new_phi =  local_machine_coords[1];
				std::vector<double> new_machine_coordinates;
				new_machine_coordinates.push_back(new_vertical_coord);
				new_machine_coordinates.push_back(new_phi);
				p_new_parent_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(r_pair.first, new_machine_coordinates));
			}
	}

	CellPropertyCollection daughter_collection = pNewCellTemp->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
	daughter_collection.RemoveProperty<TypeSixMachineProperty>();
	collection.RemoveProperty<TypeSixMachineProperty>();

	//boost::shared_ptr<TypeSixMachineProperty> p_daughter_property = boost::static_pointer_cast<TypeSixMachineProperty>(daughter_collection.GetProperty());


	collection.AddProperty(p_new_parent_property);
	daughter_collection.AddProperty(p_new_daughter_property);

	//p_parent_property=p_new_parent_property;
	//p_daughter_property=p_new_daughter_property;
	return pNewCellTemp;

}

template<unsigned DIM>
c_vector<double, DIM> NodeBasedCellPopulationWithCapsules<DIM>::GetMachineCoords(unsigned node_index,std::vector<double> machine_coordinates,c_vector<double,DIM> cell_centre,double L)
{


    Node<DIM>* p_node = this->GetNode(node_index);
    double cell_angle_theta = p_node->rGetNodeAttributes()[NA_THETA];
    double R = p_node->rGetNodeAttributes()[NA_RADIUS];


    double vertical_coordinate = machine_coordinates[0];
    double azimuthal_coordinate = machine_coordinates[1];

    double x_in_cell_frame = vertical_coordinate;


    double rho;

    if (fabs(vertical_coordinate)<L/2.0 )
    {
    	rho = R;
    }
    else if (vertical_coordinate>L/2.0)
    {
    	double z_diff=vertical_coordinate-L/2.0;
    	rho =sqrt(R*R-z_diff*z_diff);
    }
    else if (vertical_coordinate<L/2.0)
    {
    	double z_diff=vertical_coordinate+L/2.0;
    	rho =sqrt(R*R-z_diff*z_diff);
    }


    double y_in_cell_frame = rho*cos(azimuthal_coordinate);









    /*// compute angle that defines end of the axis and beginning of hemisphere
    double theta_c = atan2(R, 0.5*L);


    double theta =machine_coordinates[1];
    //TODO - DEal with second angle

    double y_in_cell_frame = DOUBLE_UNSET;
    //if (fabs(theta - 0.5*M_PI) < 0.5*M_PI - theta_c) // machine lies on upper body axis
    if (theta > theta_c && theta < M_PI-theta_c) // machine lies on upper body axis
    {
        x_in_cell_frame = R/tan(theta);
        y_in_cell_frame = R;
    }
    //else if (fabs(theta + 0.5*M_PI) < 0.5*M_PI - theta_c) // machine lies on lower body axis
    else if (theta < -theta_c && theta > -M_PI+theta_c) // machine lies on upper body axis
    {
        x_in_cell_frame = -R/tan(theta);
        y_in_cell_frame = -R;
    }
    else if (fabs(theta) < theta_c) // machine lies on right hemisphere
    {
        x_in_cell_frame = (L+sqrt(L*L-(L*L-4*R*R)*(1+tan(theta)*tan(theta))))/(2*(1+tan(theta)*tan(theta)));
        y_in_cell_frame = x_in_cell_frame*tan(theta);
    }
    else
    {
        x_in_cell_frame = (-L-sqrt(L*L-(L*L-4*R*R)*(1+tan(theta)*tan(theta))))/(2*(1+tan(theta)*tan(theta)));
        y_in_cell_frame = x_in_cell_frame*tan(theta);
    }


 */
    // Compute and store the coordinates of this machine
    c_vector<double, DIM> machine_coords = zero_vector<double>(DIM);

    if (DIM ==2)
    {
    	c_vector<double, DIM> cartesian_coords_in_cell_frame;
		cartesian_coords_in_cell_frame(0)=x_in_cell_frame;
		cartesian_coords_in_cell_frame(1)=y_in_cell_frame;

		c_matrix<double,2,2> rotation_matrix_theta = identity_matrix<double>(2);
		rotation_matrix_theta(0,0)=cos(cell_angle_theta);
		rotation_matrix_theta(0,1)=-sin(cell_angle_theta);
		rotation_matrix_theta(1,0)=sin(cell_angle_theta);
		rotation_matrix_theta(1,1)=cos(cell_angle_theta);


		machine_coords=cell_centre + prod(rotation_matrix_theta, cartesian_coords_in_cell_frame);
		//PRINT_VECTOR(machine_coords);

//    			machine_coords=cell_centre+
//
//        machine_coords[0] = cell_centre[0] + x_in_cell_frame*cos(cell_angle) - y_in_cell_frame*sin(cell_angle);
//        machine_coords[1] = cell_centre[1] + x_in_cell_frame*sin(cell_angle) + y_in_cell_frame*cos(cell_angle);
    }
    else if (DIM==3)
    {
        double cell_angle_phi = p_node->rGetNodeAttributes()[NA_PHI];
    	double z_in_cell_frame = rho*sin(azimuthal_coordinate);


    	c_vector<double, DIM> cartesian_coords_in_cell_frame;
    	cartesian_coords_in_cell_frame(0)=x_in_cell_frame;
    	cartesian_coords_in_cell_frame(1)=y_in_cell_frame;
    	cartesian_coords_in_cell_frame(2)=z_in_cell_frame;

		c_matrix<double,3,3> rotation_matrix_theta = identity_matrix<double>(3);
		rotation_matrix_theta(0,0)=cos(cell_angle_theta);
		rotation_matrix_theta(0,1)=-sin(cell_angle_theta);
		rotation_matrix_theta(1,0)=sin(cell_angle_theta);
		rotation_matrix_theta(1,1)=cos(cell_angle_theta);

		c_matrix<double,3,3> rotation_matrix_phi = identity_matrix<double>(3);
		rotation_matrix_phi(1,1)=cos(cell_angle_phi);
		rotation_matrix_phi(1,2)=-sin(cell_angle_phi);
		rotation_matrix_phi(2,1)=sin(cell_angle_phi);
		rotation_matrix_phi(2,2)=cos(cell_angle_phi);

		machine_coords=cell_centre + prod(prod(rotation_matrix_theta,rotation_matrix_phi), cartesian_coords_in_cell_frame);
    }


    return machine_coords;

}

template<unsigned DIM>
std::vector<unsigned> NodeBasedCellPopulationWithCapsules<DIM>::GetMachineData(CellPtr pCell)
{



	unsigned totalNumberMachines=0u;
	unsigned totalNumTypeL=0u;
	unsigned totalNumTypeB=0u;
	unsigned totalNumTypeH=0u;


	// Get this cell's type six machine property data
	CellPropertyCollection collection = pCell->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
	if (collection.GetSize() != 1)
	{
		EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
	}
	boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
	std::vector<std::pair<unsigned, std::vector<double>> >& r_data = p_property->rGetMachineData();

	 unsigned numMachineFiresInThisTimeStep=p_property->GetNumMachineFiresInThisTimeStep();



	// loop over machines
	for (auto& r_pair : r_data)
	{
		if (r_pair.first==1)
		{
			totalNumTypeL++;
		}
		else if (r_pair.first==2)
		{
			totalNumTypeL++;
			totalNumTypeB++;
		}
		else if (r_pair.first==3)
		{
			totalNumTypeL++;
			totalNumTypeB++;
			totalNumTypeH++;
		}
	  }



	totalNumberMachines+=r_data.size();


    std::vector<unsigned> machine_data(4);
	machine_data[0]=totalNumberMachines;
	machine_data[1]=totalNumTypeL;
	machine_data[2]=totalNumTypeB;
	machine_data[3]=totalNumTypeH;
	machine_data[4]=numMachineFiresInThisTimeStep;



    return machine_data;
}




template<unsigned DIM>
void NodeBasedCellPopulationWithCapsules<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Currently no specific parameters to output all come from parent classes

    // Call method on direct parent class
    NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

// Explicit instantiation
template class NodeBasedCellPopulationWithCapsules<1>;
template class NodeBasedCellPopulationWithCapsules<2>;
template class NodeBasedCellPopulationWithCapsules<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithCapsules)
