
#include "TypeSixMachineModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "TypeSixMachineProperty.hpp"
#include "Exception.hpp"
#include "VtkMeshWriter.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "OutputFileHandler.hpp"
#include "NodesOnlyMesh.hpp"
#include "UniformCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"
#include "NodeBasedCellPopulationWithCapsules.hpp"


template<unsigned DIM>
TypeSixMachineModifier<DIM>::TypeSixMachineModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mOutputDirectory(""),
	  mk_1(0.4),
	  mk_2(0.0),
	  mk_3(1.1),
	  mk_4(0.0),
	  mk_5(1.1),
	  mk_6(0.0),
	  mk_7(1.1),
	  mStateFire(3u)
{
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::SetMachineParametersFromGercEtAl()
{
	// Gerc et al, Cell Reports Gerc et al., 2015, Cell Reports 	12 	, 2131–2142	 2015
	// http://dx.doi.org/10.1016/j.celrep.2015.08.05

	//mk_7=1.0/30.0;
    //mk_4=0.69/2.0-mk_7;

	mk_1= 0.0005*100.0;
	mk_2=0.1;
	mk_3=0.18;
	mk_4=0.21;
	mk_5=0.1309;
	mk_6=0.10;
	mk_7=0.12;

//    mk_3=(mk_4+mk_7)/3.0;
//    mk_6=0.02*mk_7; //1.0/30.0;
//    mk_5=mk_6+mk_7;
//	mk_2=mk_6; //0.05;


	mk_1=mk_1*60.0;
    mk_2=mk_2*60.0;
    mk_3=mk_3*60.0;
    mk_4=mk_4*60.0;
    mk_5=mk_5*60.0;
    mk_6=mk_6*60.0;
    mk_7=mk_7*60.0;
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::SetContactDependentFiring()
{
	// Gerc et al, Cell Reports Gerc et al., 2015, Cell Reports 	12 	, 2131–2142	 2015
	// http://dx.doi.org/10.1016/j.celrep.2015.08.05

	mk_7=0.0;
}

template<unsigned DIM>
TypeSixMachineModifier<DIM>::~TypeSixMachineModifier()
{
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::SetOutputDirectory(std::string directory)
{
    mOutputDirectory = directory;
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::WriteVtk(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
#ifdef CHASTE_VTK



    // Store the present time as a string
    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    // Create mesh writer for VTK output
    VtkMeshWriter<DIM, DIM> mesh_writer(mOutputDirectory,
                                        "machine_results_"+time.str(),
                                        false);

    std::vector<Node<DIM>*> machine_nodes;

    // Create vector to store VTK cell data
    std::vector<double> vtk_machine_data;

    // Iterate over cell population
    unsigned machine_index = 0;


    NodeBasedCellPopulationWithCapsules<DIM>& rcapsule_pop=(static_cast<NodeBasedCellPopulationWithCapsules<DIM>&>(rCellPopulation));

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    { 
        // Get this cell's type six machine property data
        CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
        if (collection.GetSize() != 1)
        {
            EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
        }
        boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
        std::vector<std::pair<unsigned, std::vector<double>> >& r_data = p_property->rGetMachineData();
        
		Node<DIM>* p_node = rcapsule_pop.GetNodeCorrespondingToCell(*cell_iter);

		double L = p_node->rGetNodeAttributes()[NA_LENGTH];
        c_vector<double, DIM> cell_centre = p_node->rGetLocation();
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // Populate the vector of VTK data
        for (auto& r_pair : r_data)
        {
            // Populate data for VTK
            vtk_machine_data.emplace_back(r_pair.first);

            // Store the location of this machine
            c_vector<double, DIM> machine_coords=rcapsule_pop.GetMachineCoords(node_index,r_pair.second,cell_centre,L);
            machine_nodes.push_back(new Node<DIM>(machine_index, machine_coords, false));
            machine_index++;
        }
    }

    mesh_writer.AddPointData("machines", vtk_machine_data); 
    /*
     * At present, the VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK, then visualized
     * as glyphs in Paraview.
     */
    NodesOnlyMesh<DIM> temp_mesh;
    ///\todo Consider removing hardcoding of "1.0" below
    temp_mesh.ConstructNodesWithoutMesh(machine_nodes, 1.0);
    mesh_writer.WriteFilesUsingMesh(temp_mesh);

    *mpVtkMetaFile << "        <DataSet timestep=\"";
    *mpVtkMetaFile << num_timesteps;
    *mpVtkMetaFile << "\" group=\"\" part=\"0\" file=\"machine_results_";
    *mpVtkMetaFile << num_timesteps;
    *mpVtkMetaFile << ".vtu\"/>\n";

    // Tidy up
    for (unsigned i=0; i<machine_nodes.size(); i++)
    {
        delete machine_nodes[i];
    }


#endif //CHASTE_VTK
}


template<unsigned DIM>
void TypeSixMachineModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (mOutputDirectory == "")
    {
       EXCEPTION("SetOutputDirectory() must be called on a TypeSixMachineModifier before it is passed to a simulation");
    }
    
    WriteVtk(rCellPopulation);
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
#ifdef CHASTE_VTK
    // Create output files for the visualizer
    double time_now = SimulationTime::Instance()->GetTime();
    std::ostringstream time_string;
    time_string << time_now;

    if (mOutputDirectory == "")
    {
        EXCEPTION("SetOutputDirectory() must be called on TypeSixMachineModifier");
    }
    mOutputDirectory += "/machine_results_from_time_" + time_string.str();

    OutputFileHandler output_file_handler(mOutputDirectory, false);
    mpVtkMetaFile = output_file_handler.OpenOutputFile("machine_results.pvd");
    *mpVtkMetaFile << "<?xml version=\"1.0\"?>\n";
    *mpVtkMetaFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *mpVtkMetaFile << "    <Collection>\n";
#endif //CHASTE_VTK

    WriteVtk(rCellPopulation);

    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::Setk_1(double k_1)
{

	mk_1= k_1;
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::Setk_5(double k_5)
{

	mk_5= k_5;
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::Setk_2(double k_2)
{

	mk_2= k_2;
}

template<unsigned DIM>
unsigned TypeSixMachineModifier<DIM>::GetTotalNumberOfMachines(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{



	unsigned totalNumberMachines=0u;
    ///\todo Make sure the cell population is updated?
    //rCellPopulation.Update();

    // Iterate over cell population
	//PRINT_VARIABLE(rCellPopulation.GetNumRealCells());
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get this cell's type six machine property data
        CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
        if (collection.GetSize() != 1)
        {
            EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
        }
        boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
        std::vector<std::pair<unsigned, std::vector<double>> >& r_data = p_property->rGetMachineData();



        totalNumberMachines+=r_data.size();

    }


    return totalNumberMachines;
}



template<unsigned DIM>
void TypeSixMachineModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{



    ///\todo Make sure the cell population is updated?
    //rCellPopulation.Update();
    double dt = SimulationTime::Instance()->GetTimeStep();


    // Iterate over cell population and update machines
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    { 
        // Get this cell's type six machine property data
        CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
        if (collection.GetSize() != 1)
        {
        	//continue;
            EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
        }
        boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
        std::vector<std::pair<unsigned, std::vector<double>> >& r_data = p_property->rGetMachineData();
        //unsigned& r_NumMachineFiresInThisTimeStep = p_property->rGetNumMachineFiresInThisTimeStep();


        unsigned numMachineFiresInThisTimeStep=0;

        assert((mk_1 + mk_2 + mk_3 + mk_4 + mk_5 + mk_6 + mk_7 )*dt <= 1.0);

		// update existing machines

        for (auto& r_pair : r_data)
        {
	        unsigned old_state = r_pair.first;
	        unsigned new_state = old_state;
		    double r = RandomNumberGenerator::Instance()->ranf();
		
		    switch (old_state)
		    {
		        case 1u:
		            if (r < mk_2*dt)
		            {
		                new_state = 0u;
		            }
		            else if (r < (mk_2 + mk_3)*dt)
		            {
		                new_state = 2u;
		            }
		            break;
		        case 2u:
		            if (r < mk_4*dt)
		            {
		                new_state = 1u;
		            }
		            else if (r < (mk_4 + mk_5)*dt)
		            {
		                new_state = 3u;
		            }
		            break;
		        case 3u:
		            if (r < mk_6*dt)
		            {
		                new_state = 2u;
		            }
		            else if (r < (mk_6 + mk_7)*dt) // Aggressive type VI fires without any neighbour contact
		            {
		                new_state = 0u;
		                numMachineFiresInThisTimeStep++;
		            }
		            break;
//		        case 4u:
//		            if (r < mk_8*dt)
//		            {
//		                new_state = 3u;
//		            }
//		            break;
		    }
		    r_pair.first = new_state;
        }
        p_property->SetNumMachineFiresInThisTimeStep(numMachineFiresInThisTimeStep);


        //r_NumMachineFiresInThisTimeStep=
    }


		
    NodeBasedCellPopulationWithCapsules<DIM>& rcapsule_pop=(static_cast<NodeBasedCellPopulationWithCapsules<DIM>&>(rCellPopulation));


        // Iterate over cell population and create new machines randomly
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
		{
			// Get this cell's type six machine property data
			CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
			if (collection.GetSize() != 1)
			{
				EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
			}
			boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
			std::vector<std::pair<unsigned, std::vector<double>> >& r_data = p_property->rGetMachineData();


			// Create a machine?
			double r = RandomNumberGenerator::Instance()->ranf();




			if (r < mk_1*dt)
			{
                std::vector<double> machine_coordinates;

                Node<DIM>* p_node = rcapsule_pop.GetNodeCorrespondingToCell(*cell_iter);
                double L = p_node->rGetNodeAttributes()[NA_LENGTH];
                double radius = p_node->rGetNodeAttributes()[NA_RADIUS];


                double 	vertical_coordinate=(L+2.0*radius)*(RandomNumberGenerator::Instance()->ranf()-0.5);
                machine_coordinates.push_back(vertical_coordinate);

			    if (DIM==2)
			    {

			    	double r2 =RandomNumberGenerator::Instance()->ranf();

			    	double azimuthal_coordinate=M_PI;
			    	if (r2>0.5)
			    	{
			    		azimuthal_coordinate=-M_PI;
			    	}
			    	machine_coordinates.push_back(azimuthal_coordinate);
			    }
			    if (DIM ==3)
				{
			        double azimuthal_coordinate =  2*M_PI*RandomNumberGenerator::Instance()->ranf();
			        machine_coordinates.push_back(azimuthal_coordinate);

				}

				r_data.emplace_back(std::pair<unsigned, std::vector<double> >(1u, machine_coordinates));
			}


		}




    // Iterate over cells and remove machines in State 0
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            // Get this cell's type six machine property data
            CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
            if (collection.GetSize() != 1)
            {
                EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
            }
            boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
            std::vector<std::pair<unsigned, std::vector<double>> >& r_data = p_property->rGetMachineData();


            // Create a new vector to store all pairs less any we might throw away
            std::vector<std::pair<unsigned, std::vector<double>> > new_data ;
    		//new_data.reserve(r_data.size() + 1);

            for (auto& r_pair : r_data)
            {
    	        unsigned current_state = r_pair.first;

    	        if (current_state!=0u) // discard any machines in state 0
    	        {
    		        new_data.emplace_back(std::pair<unsigned, std::vector<double>>(r_pair));
    	        }

            }
    		r_data = new_data;
        }


}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
#ifdef CHASTE_VTK
    *mpVtkMetaFile << "    </Collection>\n";
    *mpVtkMetaFile << "</VTKFile>\n";
    mpVtkMetaFile->close();
#endif //CHASTE_VTK
}
    
template<unsigned DIM>
void TypeSixMachineModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class TypeSixMachineModifier<1>;
template class TypeSixMachineModifier<2>;
template class TypeSixMachineModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TypeSixMachineModifier)

