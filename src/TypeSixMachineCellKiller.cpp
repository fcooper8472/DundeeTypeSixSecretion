
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point.hpp>

#include "TypeSixMachineCellKiller.hpp"
#include "TypeSixMachineProperty.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
TypeSixMachineCellKiller<DIM>::TypeSixMachineCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
}

template<unsigned DIM>
void TypeSixMachineCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    assert (DIM == 2);
    assert(bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation)));
    NodeBasedCellPopulation<DIM>* p_population = static_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation);

    using geom_point = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;
    using geom_segment = boost::geometry::model::segment<geom_point>;

    // Iterate over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Get this cell's type six machine property data
        CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
        if (collection.GetSize() != 1)
        {
            EXCEPTION("TypeSixMachineCellKiller cannot be used unless each cell has a TypeSixMachineProperty");
        }
        boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
        std::vector<std::pair<unsigned, double> >& r_data = p_property->rGetMachineData();

		// Create a new vector to store all pairs less any we might throw away
		std::vector<std::pair<unsigned, double> > new_data;
		new_data.reserve(r_data.size() + 1);

        // Iterate over machines in this cell
        for (auto& r_pair : r_data)
        {
            // If this machine is ready to kill a cell...
            unsigned state = r_pair.first;
            if (state == 4)
            {
                // ...check if any neighbouring cells are close enough to kill...
                unsigned node_index = p_population->GetLocationIndexUsingCell(*cell_iter);
                Node<DIM>* p_node = p_population->GetNode(node_index);
                auto cell_centre = p_node->rGetLocation();
                double cell_angle = p_node->rGetNodeAttributes()[NA_ANGLE];
                double L = p_node->rGetNodeAttributes()[NA_LENGTH];
                double R = p_node->rGetNodeAttributes()[NA_RADIUS];

                double theta = r_pair.second;
                
                // We assume that the machine's angl
                double theta_c = atan2(R, 0.5*L);
                if (theta_c < 0)
                {
                    theta_c += 2*M_PI;
                }
                
                double x_in_cell_frame = DOUBLE_UNSET;
                double y_in_cell_frame = DOUBLE_UNSET;
                if (fabs(theta - 0.5*M_PI) < 0.5*M_PI - theta_c)
                {
                    x_in_cell_frame = R/tan(theta);
                    y_in_cell_frame = R; 
                }
                else if (fabs(theta - 1.5*M_PI) < 0.5*M_PI - theta_c)
                {
                    x_in_cell_frame = R/tan(theta);
                    y_in_cell_frame = -R;
                }
                else if ((theta > 2*M_PI - theta_c) || (theta < theta_c))
                {
                    x_in_cell_frame = (L+sqrt(L*L-(L*L-4*R*R)*(1+tan(theta)*tan(theta))))/(2*(1+tan(theta)*tan(theta)));
                    y_in_cell_frame = x_in_cell_frame*tan(theta); 
                }
                else  
                {
                	x_in_cell_frame = (-L-sqrt(L*L-(L*L-4*R*R)*(1+tan(theta)*tan(theta))))/(2*(1+tan(theta)*tan(theta)));
                    y_in_cell_frame = x_in_cell_frame*tan(theta); 
                }

	            // Compute and store the coordinates of this machine            
	            c_vector<double, DIM> machine_coords;
	            machine_coords[0] = cell_centre[0] + x_in_cell_frame*cos(cell_angle) - y_in_cell_frame*sin(cell_angle);
	            machine_coords[1] = cell_centre[1] + x_in_cell_frame*sin(cell_angle) + y_in_cell_frame*cos(cell_angle);
                
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
                    auto neighbour_location = p_neighbour->rGetLocation();
				    const double neighbour_angle = p_neighbour->rGetNodeAttributes()[NA_ANGLE];
				    const double neighbour_length = p_neighbour->rGetNodeAttributes()[NA_LENGTH];
				
				    geom_point machine_point(machine_coords[0], machine_coords[1]);
								
				    geom_point neighbour_end_1(neighbour_location[0] + 0.5 * neighbour_length * cos(neighbour_angle),
				                               neighbour_location[1] + 0.5 * neighbour_length * sin(neighbour_angle));
				
				    geom_point neighbour_end_2(neighbour_location[0] - 0.5 * neighbour_length * cos(neighbour_angle),
				                               neighbour_location[1] - 0.5 * neighbour_length * sin(neighbour_angle));
				
				    geom_segment neighbour_axis(neighbour_end_1, neighbour_end_2);
				
				    double distance_to_neighbour = boost::geometry::distance(machine_point, neighbour_axis);
                    
                    if (distance_to_neighbour > R)
                    {
		                new_data.emplace_back(std::pair<unsigned, double>(r_pair));
		            }
		            else
		            {      
                        // Kill this neighbouring cell
                        p_population->GetCellUsingLocationIndex(*it)->StartApoptosis();
                        break;
                        
	                    // Note: In this case we don't store this machine's data in new_data,
	                    // since the machine is assumed to be destroyed upon killing the 
	                    // neighbouring cell.
                    }
                }
            }
            else
            {
		        new_data.emplace_back(std::pair<unsigned, double>(r_pair));
		    }
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
