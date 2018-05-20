
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point.hpp>

#include "TypeSixMachineCellKiller.hpp"
#include "TypeSixMachineProperty.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"

#include "NodeBasedCellPopulationWithCapsules.hpp"

template<unsigned DIM>
TypeSixMachineCellKiller<DIM>::TypeSixMachineCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
}

template<unsigned DIM>
void TypeSixMachineCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    assert (DIM >1);
    assert(bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation)));
    NodeBasedCellPopulation<DIM>* p_population = static_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation);

    using geom_point = boost::geometry::model::point<double, DIM, boost::geometry::cs::cartesian>;
    using geom_segment = boost::geometry::model::segment<geom_point>;

    // Iterate over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = p_population->Begin();
         cell_iter != p_population->End();
         ++cell_iter)
    {
        // Get this cell's type six machine property data
        CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
        if (collection.GetSize() != 1)
        {
            EXCEPTION("TypeSixMachineCellKiller cannot be used unless each cell has a TypeSixMachineProperty");
        }
        boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
        std::vector<std::pair<unsigned, std::vector<double>> >& r_data = p_property->rGetMachineData();

		// Create a new vector to store all pairs less any we might throw away
        std::vector<std::pair<unsigned, std::vector<double>> > new_data ;

        // Iterate over machines in this cell
        for (auto& r_pair : r_data)
        {
            // If this machine is ready to kill a cell...
            unsigned state = r_pair.first;
            if (state == 3u)
            {
                // ...check if any neighbouring cells are close enough to kill...
                unsigned node_index = p_population->GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, DIM> cell_centre = p_population->GetNodeCorrespondingToCell(*cell_iter)->rGetLocation();

	            NodeBasedCellPopulationWithCapsules<DIM>* p_capsule_pop=(dynamic_cast<NodeBasedCellPopulationWithCapsules<DIM>*>(p_population));
	            Node<DIM>* p_node = p_capsule_pop->GetNodeCorrespondingToCell(*cell_iter);
	           	double L = p_node->rGetNodeAttributes()[NA_LENGTH];
	            c_vector<double, DIM> machine_coords=p_capsule_pop->GetMachineCoords(node_index,r_pair.second,cell_centre,L);

                // Store the node indices corresponding to neighbouring cells
                double neighbourhood_radius = 3.0*L;
                std::set<unsigned> neighbours = p_population->GetNodesWithinNeighbourhoodRadius(node_index, neighbourhood_radius);

                // Iterate over neighbouring cells
                for (std::set<unsigned>::iterator it = neighbours.begin();
                     it != neighbours.end();
                     ++it)
                {
                    Node<DIM>* p_neighbour = p_population->GetNode(*it);
                    
                    // Compute distance between (X,Y) and this neighbouring cell's line segment
				    const double neighbour_angle = p_neighbour->rGetNodeAttributes()[NA_THETA];
				    const double neighbour_length = p_neighbour->rGetNodeAttributes()[NA_LENGTH];
				    double R = p_neighbour->rGetNodeAttributes()[NA_RADIUS];

				    double distance_to_neighbour =0.0;
				    if (DIM==2u)
				    {
				    	geom_point machine_point(machine_coords[0], machine_coords[1]);

				    	c_vector<double, DIM> neighbour_plus_end_local;
				    	neighbour_plus_end_local(0)=neighbour_length/2.0;
				    	neighbour_plus_end_local(1)=0;
				    	c_vector<double, DIM> neighbour_min_end_local;
				    	neighbour_min_end_local(0)=-neighbour_length/2.0;
				    	neighbour_min_end_local(1)=0;

						c_matrix<double,2,2> rotation_matrix_theta = identity_matrix<double>(2);
						rotation_matrix_theta(0,0)=cos(neighbour_angle);
						rotation_matrix_theta(0,1)=-sin(neighbour_angle);
						rotation_matrix_theta(1,0)=sin(neighbour_angle);
						rotation_matrix_theta(1,1)=cos(neighbour_angle);

						c_vector<double, DIM> neighbour_plus_end_global=p_neighbour->rGetLocation() + prod(rotation_matrix_theta, neighbour_plus_end_local);
						c_vector<double, DIM> neighbour_min_end_global=p_neighbour->rGetLocation() + prod(rotation_matrix_theta, neighbour_min_end_local);

				    	geom_point neighbour_end_1(neighbour_plus_end_global[0],neighbour_plus_end_global[1]);
				    	geom_point neighbour_end_2(neighbour_min_end_global[0],neighbour_min_end_global[1]);
				    	geom_segment neighbour_axis(neighbour_end_1, neighbour_end_2);
				    	distance_to_neighbour = boost::geometry::distance(machine_point, neighbour_axis);
				    }
				    else if (DIM==3u)
				    {
					    const double neighbour_phi = p_neighbour->rGetNodeAttributes()[NA_PHI];

				    	geom_point machine_point(machine_coords[0], machine_coords[1], machine_coords[2]);

				    	c_vector<double, DIM> neighbour_plus_end_local;
						neighbour_plus_end_local(0)=neighbour_length/2.0;
						neighbour_plus_end_local(1)=0.0;
						neighbour_plus_end_local(2)=0.0;

						c_vector<double, DIM> neighbour_min_end_local;
						neighbour_min_end_local(0)=-neighbour_length/2.0;
						neighbour_min_end_local(1)=0.0;
						neighbour_min_end_local(2)=0.0;

						c_matrix<double,3,3> rotation_matrix_theta = identity_matrix<double>(3);
						rotation_matrix_theta(0,0)=cos(neighbour_angle);
						rotation_matrix_theta(0,1)=-sin(neighbour_angle);
						rotation_matrix_theta(1,0)=sin(neighbour_angle);
						rotation_matrix_theta(1,1)=cos(neighbour_angle);

						c_matrix<double,3,3> rotation_matrix_phi = identity_matrix<double>(3);
						rotation_matrix_phi(1,1)=cos(neighbour_phi);
						rotation_matrix_phi(1,2)=-sin(neighbour_phi);
						rotation_matrix_phi(2,1)=sin(neighbour_phi);
						rotation_matrix_phi(2,2)=cos(neighbour_phi);

						c_vector<double, DIM> neighbour_plus_end_global=p_neighbour->rGetLocation() + prod(prod(rotation_matrix_theta,rotation_matrix_phi), neighbour_plus_end_local);

						c_vector<double, DIM> neighbour_min_end_global=p_neighbour->rGetLocation() + prod(prod(rotation_matrix_theta,rotation_matrix_phi), neighbour_min_end_local);

						geom_point neighbour_end_1(neighbour_plus_end_global[0],neighbour_plus_end_global[1],neighbour_plus_end_global[2]);
						geom_point neighbour_end_2(neighbour_min_end_global[0],neighbour_min_end_global[1], neighbour_min_end_global[2]);

						geom_segment neighbour_axis(neighbour_end_1, neighbour_end_2);

						distance_to_neighbour = boost::geometry::distance(machine_point, neighbour_axis);
				    }
					else
					{
				            EXCEPTION("TypeSixMachineCellKiller cannot be used unless each cell has a TypeSixMachineProperty");
					}
                    
                    if (distance_to_neighbour > R || p_population->GetCellUsingLocationIndex(*it)->HasApoptosisBegun())
                    {
		                //new_data.emplace_back(std::pair<unsigned, double>(r_pair));
		            }
		            else
		            {
		            	TRACE("StartingApoptosis");
                        // Kill this neighbouring cell
		                r_pair.first=0u;

                        p_population->GetCellUsingLocationIndex(*it)->StartApoptosis();
                        // Note: In this case we don't store this machine's data in new_data,
	                    // since the machine is assumed to be destroyed upon killing the 
	                    // neighbouring cell.
                    }
                }
            }
            else
            {
		    }
            new_data.emplace_back(std::pair<unsigned, std::vector<double>>(r_pair));
        }
		r_data = new_data;
    }
}

template<unsigned DIM>
void TypeSixMachineCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class TypeSixMachineCellKiller<1>;
template class TypeSixMachineCellKiller<2>;
template class TypeSixMachineCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TypeSixMachineCellKiller)
