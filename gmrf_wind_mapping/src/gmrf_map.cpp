#include "gmrf_wind_mapping/gmrf_map.h"
#include <cstdio>   // Necesario para fprintf
#include <iostream>
#include <iomanip>  // Necesario para std::setprecision


/*---------------------------------------------------------------
                        Constructor
  ---------------------------------------------------------------*/
CGMRF_map::CGMRF_map(const TOccupancyMap& oc_map, 
                     float cell_size, 
                     double m_lambdaPrior_reg,
                     double m_lambdaPrior_flux_conservation, 
                     double m_lambdaPrior_obstacles,
                     double m_lambdaObservations,
                     bool verbose,
                     bool estimateTiming=false)
{
    // Set Verbose level
    this->verbose = verbose;
    this->estimateTiming = estimateTiming;

    try
    {
        // Copy params to internal variables
        m_Ocgridmap = oc_map;     // Occupancy gridMap ( from ROS2 MapServer or other sources)
        m_resolution = cell_size; // Desired resolution to build the GMRF (m)
        lambdaPrior_reg = m_lambdaPrior_reg;
        lambdaPrior_flux_conservation = m_lambdaPrior_flux_conservation;
        lambdaPrior_obstacles = m_lambdaPrior_obstacles;
        lambdaObservations = m_lambdaObservations;
        // Compute SQRT values
        //lambdaPrior_reg_sqrt = std::sqrt(lambdaPrior_reg);
        //lambdaPrior_flux_conservation_sqrt = std::sqrt(lambdaPrior_flux_conservation);
        //lambdaPrior_obstacles_sqrt = std::sqrt(lambdaPrior_obstacles);

        // Set initial GMRF dimensions as the OccupancyMap (in meters)
        double x_min = oc_map.origin_x;
        double x_max = oc_map.origin_x + oc_map.width * oc_map.resolution;
        double y_min = oc_map.origin_y;
        double y_max = oc_map.origin_y + oc_map.height * oc_map.resolution;

        // Adjust size to complaint with the desired resolution (m_resolution):
        m_x_min = m_resolution * round(x_min / m_resolution);
        m_y_min = m_resolution * round(y_min / m_resolution);
        m_x_max = m_resolution * round(x_max / m_resolution);
        m_y_max = m_resolution * round(y_max / m_resolution);

        // Now the number of cells should be integers:
        m_size_x = round((m_x_max - m_x_min) / m_resolution);
        m_size_y = round((m_y_max - m_y_min) / m_resolution);
        N = m_size_x * m_size_y;    // Number of cells in the GMRF


        // 1. INIT RANDOM FIELD AND CREATE CONNEXIONS BETWEEN NODES
        //-----------------------------------------------------------
        std::cerr << "[GMRF_MAP] Generating GMRF for 2D WIND estimation..." << std::endl;

        // 1. Init the map container (2N cells)
        //-------------------------
        TRandomFieldCell init_cell;
        init_cell.mean = 0.0;
        init_cell.std = 0.0;
        m_map.assign(2 * N, init_cell); // Since we have Wx and Wy, we refer to them as: Wx in the range [0,N-1], Wy in the range [N,2N-1]

        if (verbose)
        {
            std::cerr << "--------------------------------" << std::endl;
            std::cerr << "[CGMRF] GMRF created:" << std::endl;
            std::cerr << "[CGMRF] Using OccupancyGrid with limits: x=(" << x_min << "," << x_max << ") [m] and y=(" << y_min << "," << y_max << std::endl;
            std::cerr << "[CGMRF] Using OccupancyGrid with cell size: (" << m_Ocgridmap.width << "," << m_Ocgridmap.height << ") cells with cell_size " << std::fixed << std::setprecision(2) << m_Ocgridmap.resolution << "[m]" << std::endl;
            std::cerr << "[CGMRF] GMRF limits: x=(" << std::fixed << std::setprecision(2) << m_x_min << "," << m_x_max << ")[m] y=(" << m_y_min << "," << m_y_max << ")[m]" << std::endl;
            std::cerr << "[CGMRF] GMRF cell size: (" << m_size_x << "," << m_size_y << ") cells with cell_size " << std::fixed << std::setprecision(2) << m_resolution << "[m]" << std::endl;
            std::cerr << "[CGMRF] GMRF with " << N << " cells and " << m_map.size() << " nodes" << std::endl;
            std::cerr <<  "--------------------------------" << std::endl;
        }

        // 2. Memory Reservation (seepUp)
        //-------------------------------
        // Approximate number of static factors
        nPriorFactors = 2 * ((m_size_x - 1) * m_size_y + m_size_x * (m_size_y - 1));
        nObsFactors = 0;
        nFactors = nPriorFactors + nObsFactors;
        if (verbose)
            std::cerr <<  "[CGMRF] Reserving Memory for Prior-factors" << std::endl;

        // Reserve memory for the Jacobian and Information matrix
        J.clear();
        J.reserve(5 * nFactors); // Each factor accounts for the connectivity of multiple nodes -->multiple entries
        Lambda.clear();
        Lambda.reserve(nFactors); // Diagonal -> 1 entry for each factor

       
        // 3. Build Prior Factors (just once if the map is static)
        //----------------------------------------------------------
        // Given the possibility of using different cell_sizes for the GMRF 
        // and the provided OccupancyGrid map, we need to employ the OccupancyGrid
        // map to determine the real interconnections between cells and then 
        // create the factors between nodes in the GMRF.
        size_t count = 0;
        for (size_t j = 0; j < N; j++) // For each cell in the GMRF
        {
            // Get cell_x and cell_y in the GMRF representation
            size_t jx, jy;
            id2cellxy(j, jx, jy);

            if (!is_cell_free(j))
            {
                // Cell is occupied in the provided OccupancyGrid --> Estimation here has no sense
                // Since we cannot remove this cell, we will force Wx=0 and Wy=0
                // Wx(j) = 0               
                Eigen::Triplet<double> J_entry(count, j, 1);
                J.push_back(J_entry);
                factor_types.push_back({count, FactorType::Obstacle});
                //Eigen::Triplet<double> lambda_entry(count, count, lambdaPrior_obstacles);
                //Lambda.push_back(lambda_entry);
                count++;
                // Wy(j) = 0
                Eigen::Triplet<double> J_entry2(count, j + N, 1);               
                J.push_back(J_entry2);
                factor_types.push_back({count, FactorType::Obstacle});
                //Eigen::Triplet<double> lambda_entry2(count, count, lambdaPrior_obstacles);
                //Lambda.push_back(lambda_entry2);
                count++;
            }

            // 3.1 Regularization and Obstacle based Factors with the right node: (j <--> j+1)
            //----------------------------------------------------------------------------------
            if (jx < (m_size_x - 1))
            {
                if (is_cell_free(j) && is_cell_free(size_t(j + 1)))
                {
                    if (check_connectivity_between2cells(j, size_t(j + 1)))
                    {
                        // Regularization Factor: cells j and j+1 should have similar wind values
                        // Wx range [1,N]
                        Eigen::Triplet<double> J_entry1(count, j, 1);
                        Eigen::Triplet<double> J_entry2(count, size_t(j + 1), -1);
                        J.push_back(J_entry1);
                        J.push_back(J_entry2);
                        factor_types.push_back({count, FactorType::Regularization});
                         //Eigen::Triplet<double> lambda_entry(count, count, lambdaPrior_reg);
                        //Lambda.push_back(lambda_entry);
                        count++;

                        // Wy range [N+1,2N]                        
                        Eigen::Triplet<double> J_entry3(count, j + N, 1);
                        Eigen::Triplet<double> J_entry4(count, size_t(j + N + 1), -1);
                        J.push_back(J_entry3);
                        J.push_back(J_entry4);
                        factor_types.push_back({count, FactorType::Regularization});
                        //Eigen::Triplet<double> lambda_entry2(count, count, lambdaPrior_reg);
                        //Lambda.push_back(lambda_entry2);
                        count++;
                    }
                    else
                    {
                        // An obstacle is in between cells j and j+1
                        // 1. Do not create a regularization link
                        // 2. Force Wx=0 at both cells
                        // Wx(j) = 0
                        Eigen::Triplet<double> J_entry(count, j, 1);
                        J.push_back(J_entry);
                        factor_types.push_back({count, FactorType::Obstacle});
                        //Eigen::Triplet<double> lambda_entry(count, count, lambdaPrior_obstacles);
                        //Lambda.push_back(lambda_entry);
                        count++;
                        
                        // Wx(j+1) = 0
                        Eigen::Triplet<double> J_entry2(count, j + 1, 1);
                        J.push_back(J_entry2);
                        factor_types.push_back({count, FactorType::Obstacle});
                        //Eigen::Triplet<double> lambda_entry2(count, count, lambdaPrior_obstacles);
                        //Lambda.push_back(lambda_entry2);
                        count++;
                    }
                }
                else if (is_cell_free(j))
                {
                    // Cell j+1 is occupied. Force Wx=0 at j
                    Eigen::Triplet<double> J_entry(count, j, 1);
                    J.push_back(J_entry);
                    factor_types.push_back({count, FactorType::Obstacle});
                    count++;
                    //Eigen::Triplet<double> lambda_entry(count, count, lambdaPrior_obstacles);
                    //Lambda.push_back(lambda_entry);
                }
                else if (is_cell_free(j + 1))
                {
                    // Cell j is occupied. Force Wx=0 at j+1
                    Eigen::Triplet<double> J_entry(count, j + 1, 1);
                    J.push_back(J_entry);
                    factor_types.push_back({count, FactorType::Obstacle});
                    count++;
                    //Eigen::Triplet<double> lambda_entry(count, count, lambdaPrior_obstacles);
                    //Lambda.push_back(lambda_entry);
                }
                // else --> Both cells occupied -> Do nothing!
            }

            // 3.2 Regularization and Obstacle based Factors with the upper node: (j <--> j+m_size_x)
            //----------------------------------------------------------------------------------------
            if (jy < (m_size_y - 1))
            {
                if (is_cell_free(j) && is_cell_free(j + m_size_x))
                {
                    if (check_connectivity_between2cells(j, j + m_size_x))
                    {
                        // Regularization Factor: cells j and j+m_size_x should have similar wind values
                        //  Wx range [1,N]
                        Eigen::Triplet<double> J_entry1(count, j, 1);
                        Eigen::Triplet<double> J_entry2(count, j + m_size_x, -1);
                        J.push_back(J_entry1);
                        J.push_back(J_entry2);
                        factor_types.push_back({count, FactorType::Regularization});
                        count++;
                        //Eigen::Triplet<double> lambda_entry(count, count, lambdaPrior_reg);
                        //Lambda.push_back(lambda_entry);

                        // Wy range [N+1,2N]
                        Eigen::Triplet<double> J_entry3(count, j + N, 1);
                        Eigen::Triplet<double> J_entry4(count, j + N + m_size_x, -1);
                        J.push_back(J_entry3);
                        J.push_back(J_entry4);
                        factor_types.push_back({count, FactorType::Regularization});
                        count++;
                        // Eigen::Triplet<double> lambda_entry2(count, count, lambdaPrior_reg);
                        // Lambda.push_back(lambda_entry2);
                    }
                    else
                    {
                        // An obstacle is in between both cells
                        // 1. Do not create a regularization link
                        // 2. Force Wy=0 at both cells
                        // Force Wy=0 at j
                        Eigen::Triplet<double> J_entry(count, j + N, 1);                        
                        J.push_back(J_entry);
                        factor_types.push_back({count, FactorType::Obstacle});
                        count++;
                        // Eigen::Triplet<double> lambda_entry(count, count, lambdaPrior_obstacles);
                        // Lambda.push_back(lambda_entry);
                        
                        // Force Wy=0 at j+m_size_x
                        Eigen::Triplet<double> J_entry2(count, j + N + m_size_x, 1);
                        J.push_back(J_entry2);
                        factor_types.push_back({count, FactorType::Obstacle});
                        count++;
                        // Eigen::Triplet<double> lambda_entry2(count, count, lambdaPrior_obstacles);
                        // Lambda.push_back(lambda_entry2);                        
                    }
                }
                else if (is_cell_free(j))
                {
                    // Cell j+m_size_x is occupied. Force Wy=0 at j
                    Eigen::Triplet<double> J_entry(count, j + N, 1);                    
                    J.push_back(J_entry);
                    factor_types.push_back({count, FactorType::Obstacle});
                    count++;
                    //Eigen::Triplet<double> lambda_entry(count, count, lambdaPrior_obstacles);
                    //Lambda.push_back(lambda_entry);
                }
                else if (is_cell_free(j + m_size_x))
                {
                    // Cell j is occupied. Force Wy=0 at j+m_size_x
                    Eigen::Triplet<double> J_entry2(count, j + N + m_size_x, 1);
                    J.push_back(J_entry2);
                    factor_types.push_back({count, FactorType::Obstacle});
                    count++;
                    // Eigen::Triplet<double> lambda_entry2(count, count, lambdaPrior_obstacles);
                    //Lambda.push_back(lambda_entry2);
                }
                // else --> Both cells occupied -> Do nothing!
            }


            // 3.3 Factors for Law of Mass Conservation (avoid borders of map)
            //----------------------------------------------------------------------------------------
            if (is_cell_free(j) && jx > 0 && jx < m_size_x - 1 && jy > 0 && jy < m_size_y - 1)
            {
                // As soon as any of its 8 clossest neighbour cells is free, set the factor
                bool set = false;
                if (is_cell_free(j - 1))
                {
                    Eigen::Triplet<double> J_entry(count, j - 1, -1);
                    J.push_back(J_entry);
                    set = true;
                }
                if (is_cell_free(j + 1))
                {
                    Eigen::Triplet<double> J_entry(count, j + 1, 1);
                    J.push_back(J_entry);
                    set = true;
                }
                if (is_cell_free(j - m_size_x))
                {
                    Eigen::Triplet<double> J_entry(count, j + N - m_size_x, -1);
                    J.push_back(J_entry);
                    set = true;
                }
                if (is_cell_free(j + m_size_x))
                {
                    Eigen::Triplet<double> J_entry(count, j + N + m_size_x, 1);
                    J.push_back(J_entry);
                    set = true;
                }

                // Diagonals
                if (is_cell_free(j + m_size_x - 1))
                {
                    Eigen::Triplet<double> J_entry(count, j + m_size_x - 1, -0.5);
                    Eigen::Triplet<double> J_entry2(count, j + m_size_x - 1 + N, 0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }
                if (is_cell_free(j + m_size_x + 1))
                {
                    Eigen::Triplet<double> J_entry(count, j + m_size_x + 1, 0.5);
                    Eigen::Triplet<double> J_entry2(count, j + m_size_x + 1 + N, 0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }
                if (is_cell_free(j - m_size_x + 1))
                {
                    Eigen::Triplet<double> J_entry(count, j - m_size_x + 1, 0.5);
                    Eigen::Triplet<double> J_entry2(count, j - m_size_x + 1 + N, -0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }
                if (is_cell_free(j - m_size_x - 1))
                {
                    Eigen::Triplet<double> J_entry(count, j - m_size_x - 1, -0.5);
                    Eigen::Triplet<double> J_entry2(count, j - m_size_x - 1 + N, -0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }

                // If factor is to be set...
                if (set)
                {
                    factor_types.push_back({count, FactorType::FluxConservation});
                    //Eigen::Triplet<double> lambda_entry(count, count, lambdaPrior_flux_conservation);
                    //Lambda.push_back(lambda_entry);
                    count++;
                    set = false;
                }
            }
        } // end for setting factors


        // Total Number of Factors
        //-------------------------
        nPriorFactors = count;
        nFactors = nPriorFactors + nObsFactors;
        activeObs.clear();
        std::cerr <<  "[CGMRF] Initialization Complete: " << nFactors << " factors for a map size of 2N=" << m_map.size() << " nodes" << std::endl;

        
        // DEBUG: Save to file
        //-------------------
        /*
            Eigen::VectorXd y_empty;
            y_empty.resize(nFactors);
            y_empty.fill(0.0);
            save_grmf_factor_graph(J,Lambda,y_empty);
            */

            // DEGUB : ADD FIXED OBSERVATION
            /*
            TobservationGMRF new_obs;
            const int cellIdx = xy2idx( 3.0, 3.0 );
            new_obs.cell_idx = cellIdx;
            new_obs.windX = 0.10;
            new_obs.windY = 1.0;
            new_obs.lambda = 13;
            new_obs.time_invariant = false;		//Default behaviour, the obs will lose weight with time.
            std::cerr <<  "[GMRF] DEMO obs: Wx = %.2f m/s Wy = %.2f m/s at cell %lu\n\n", new_obs.windX,new_obs.windY,new_obs.cell_idx);
            activeObs.push_back(new_obs);
            nObsFactors += 2;    //we add 2 factors for each observation to account for Wx and Wy components
        */
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-Constructor] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
    }
}

/*---------------------------------------------------------------
                        Destructor
  ---------------------------------------------------------------*/
CGMRF_map::~CGMRF_map()
{
}


void CGMRF_map::update_lambdas(double m_lambdaPrior_reg,
                        double m_lambdaPrior_flux_conservation, 
                        double m_lambdaPrior_obstacles,
                        double m_lambdaObservations)
{
    // Update internal variables
    lambdaPrior_reg = m_lambdaPrior_reg;
    lambdaPrior_flux_conservation = m_lambdaPrior_flux_conservation;
    lambdaPrior_obstacles = m_lambdaPrior_obstacles;
    lambdaObservations = m_lambdaObservations;
}

void CGMRF_map::read_lambdas(double &m_lambdaPrior_reg, double &m_lambdaPrior_flux_conservation, double &m_lambdaPrior_obstacles, double &m_lambdaObservations)
{
    m_lambdaPrior_reg = lambdaPrior_reg;
    m_lambdaPrior_flux_conservation = lambdaPrior_flux_conservation;
    m_lambdaPrior_obstacles = lambdaPrior_obstacles;
    m_lambdaObservations = lambdaObservations;
}

/*---------------------------------------------------------------
                        Cell index transformations
  ---------------------------------------------------------------*/
// Get cell(x,y) in the GMRF representation from the general index in the array
void CGMRF_map::id2cellxy(size_t id, size_t& cell_x, size_t& cell_y)
{
    cell_x = id % m_size_x;
    cell_y = (size_t)floor(id / m_size_x);
}

// Get pose x,y (in meters) (in the GMRF representation) from the general index in the array
void CGMRF_map::id2xy(size_t id, double& x, double& y)
{
    size_t cell_x, cell_y;
    id2cellxy(id, cell_x, cell_y);

    x = m_x_min + (cell_x * m_resolution) + (m_resolution / 2);
    y = m_y_min + (cell_y * m_resolution) + (m_resolution / 2);
}


/*---------------------------------------------------------------
             Check if a cell is free of obstacles
  ---------------------------------------------------------------*/
// Check at OccupancyMap level, if a cell is free of obstacles (checking the cell center at GMRF resolution)
bool CGMRF_map::is_cell_free(size_t id_gmrf)
{
    // The pose x,y (meters) of cell center
    double cell_1_x, cell_1_y;
    id2xy(id_gmrf, cell_1_x, cell_1_y);

    // Get corresponding cell_idx in the Occupancy Gridmap
    // IMPORTANT --> Use the resolution and size of Occupancy Gridmap (not the GMRF)
    int id_oc;
    id_oc = static_cast<int>((cell_1_x - m_Ocgridmap.origin_x) / m_Ocgridmap.resolution);                      // x component
    id_oc += static_cast<int>((cell_1_y - m_Ocgridmap.origin_y) / m_Ocgridmap.resolution) * m_Ocgridmap.width; // y component

    /* DEBUG
    double cellx,celly;
    cellx = m_Ocgridmap.info.origin.position.x + (id_oc % m_Ocgridmap.info.width)*m_Ocgridmap.info.resolution + m_Ocgridmap.info.resolution/2;
    celly = m_Ocgridmap.info.origin.position.y + ((size_t) floor(id_oc/m_Ocgridmap.info.width))*m_Ocgridmap.info.resolution +
    m_Ocgridmap.info.resolution/2;
    */

    try
    {
        // Check occupancy
        if (m_Ocgridmap.data[id_oc] >= 50.0)
        {
            // std::cerr <<  "[GMRF] OCCUPIED %lu = (%.2f,%.2f) --> %lu in Occ",idx_1_gmrf,cell_1_x,cell_1_y, idx_1_oc);
            return false;
        }
        else
        {
            // std::cerr <<  "[GMRF] FREE %lu = (%.2f,%.2f) --> %lu in Occ",idx_1_gmrf,cell_1_x,cell_1_y, idx_1_oc);
            return true;
        }
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-is_cell_free] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
        return false;
    }
}


/*---------------------------------------------------------------
             Check cell interconnectivity
  ---------------------------------------------------------------*/
// Check at OccupancyMap level, if two cells are interconnected, that is there exists no obstacles in between them.
bool CGMRF_map::check_connectivity_between2cells(size_t idx_1_gmrf, size_t idx_2_gmrf)
{
    try
    {
        // Ge poses (x,y) of the cell centers in GMRF map
        double cell_1_x, cell_1_y, cell_2_x, cell_2_y;
        id2xy(idx_1_gmrf, cell_1_x, cell_1_y);
        id2xy(idx_2_gmrf, cell_2_x, cell_2_y);

        // Get corresponding cell_idx in the Occupancy Gridmap
        // IMPORTANT --> Use the resolution and size of Occupancy Gridmap (not the GMRF)
        int idx_1_oc, idx_2_oc;
        idx_1_oc = static_cast<int>((cell_1_x - m_Ocgridmap.origin_x) / m_Ocgridmap.resolution);                        // x component
        idx_1_oc += static_cast<int>((cell_1_y - m_Ocgridmap.origin_y) / m_Ocgridmap.resolution) * m_Ocgridmap.width;   // y component
        idx_2_oc = static_cast<int>((cell_2_x - m_Ocgridmap.origin_x) / m_Ocgridmap.resolution);                        // x component
        idx_2_oc += static_cast<int>((cell_2_y - m_Ocgridmap.origin_y) / m_Ocgridmap.resolution) * m_Ocgridmap.width;   // y component

        // check if cells are in the same row of the GMRF map
        bool horizontal = false;
        if (idx_2_gmrf == idx_1_gmrf + 1)
            horizontal = true;

        // Check that a straigh line between both cells centers is free of obstacles
        bool connected = true;
        if (horizontal)
        {
            for (size_t p = idx_1_oc; p < idx_2_oc; p++)
            {
                if (m_Ocgridmap.data[p] >= 60.0)
                {
                    connected = false;
                    break;
                }
            }
        }
        else
        {
            for (size_t p = idx_1_oc; p < idx_2_oc; p += m_Ocgridmap.width)
            {
                if (m_Ocgridmap.data[p] >= 50.0)
                {
                    connected = false;
                    break;
                }
            }
        }
        return connected;
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-check_connectivity_between2cells] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
        return false;
    }
}


/*---------------------------------------------------------------
             Insert New Wind Observation
---------------------------------------------------------------*/
void CGMRF_map::insertObservation_GMRF(double wind_speed, double wind_direction, double x_pos, double y_pos, double lambdaObs)
{
    try
    {
        auto add_obs = [this](const TobservationGMRF& observation)
        {
            if(observation.cell_idx <0 || observation.cell_idx > N)
            {
                std::cerr << "[GMRF-MAP] Observation is outside of the map!" << std::endl;
                return;
            }
            activeObs.push_back(observation);
        };

        // Get cell indexes
        const int cellIdx = xy2idx(x_pos, y_pos);

        // Fill new Observation
        // The wind vector provided is already the DownWind direction in the map reference system
        if (x_pos <= m_x_min || x_pos >= m_x_max || y_pos <= m_y_min || y_pos >= m_y_max || !is_cell_free(cellIdx))
            return;

        TobservationGMRF new_obs;
        new_obs.cell_idx = cellIdx;
        new_obs.windX = wind_speed * cos(wind_direction);
        new_obs.windY = wind_speed * sin(wind_direction);
        new_obs.lambda = lambdaObs;
        new_obs.time_invariant = true; // Default behaviour, the obs will not lose weight with time.
        if (verbose)
            std::cerr << "[GMRF-MAP] New obs: Wx = " << new_obs.windX << " m/s Wy = " << new_obs.windY << " m/s" << std::endl;

        // Add Observation to GMRF
        add_obs(new_obs);
        nObsFactors += 2; // we add 2 factors foe each observation to account for Wx and Wy components
        
        /*
            // NOTE --> Create 4 observations to expand a bit the measurement impact, replicating the content to neighbour cells
            if (is_cell_free(cellIdx - 1))
            {
                new_obs.cell_idx = cellIdx - 1;
                add_obs(new_obs);
                nObsFactors += 2;
            }
            if (is_cell_free(cellIdx - m_size_x))
            {
                new_obs.cell_idx = cellIdx - m_size_x;
                add_obs(new_obs);
                nObsFactors += 2;
            }
            if (is_cell_free(cellIdx - m_size_x - 1))
            {
                new_obs.cell_idx = cellIdx - m_size_x - 1;
                add_obs(new_obs);
                nObsFactors += 2;
            }
        */
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-insertObservation_GMRF] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
    }
}


void CGMRF_map::clearObservations_GMRF()
{
    try
    {
        activeObs.clear();
        nObsFactors = 0;
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-clearObservations_GMRF] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
    }
}


// ----------------------------------------------------------------------
// ------------ LOG MARGINAL LIKELIHOOD OPTIMIZATION LOOP ---------------
// ----------------------------------------------------------------------
void CGMRF_map::optimize_LML(bool performLMLOptimization, int maxIterations, double learningRate, double LML_threshold)
{
    if (verbose) std::cerr << "[GMRF] Starting GMRF Estimation..." << (performLMLOptimization ? " (LML Optimization Active)" : " (Fixed Hyperparameters)") << std::endl;
    
    // Always run at least one MAP estimation, even if not optimizing
    if (!performLMLOptimization) {
        MAP_estimation_GMRF();
        return; // MAP done, exit
    }

    // LML Optimization Loop -> We optimize the hyperparameters (lambdas) by maximizing the Log Marginal Likelihood
    //-----------------------------------------------------------
    double currentLML = -std::numeric_limits<double>::infinity();
    double prevLML = -std::numeric_limits<double>::infinity();
    
    // Define the lambda vector for easy update and gradient calculation
    // Order: [0]Observations, [1]Regularization, [2]Flux, [3]Obstacle
    Eigen::Vector4d lambdaVec;
    lambdaVec << lambdaObservations, lambdaPrior_reg, lambdaPrior_flux_conservation, lambdaPrior_obstacles;
    
    // Momentum and Velocity state (for optimization)
    Eigen::Vector4d velocity = Eigen::Vector4d::Zero();
    double momentum = 0.9;
    Eigen::Vector4d gradientVec;
    
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        if (verbose) std::cerr << "\n[LML Iteration " << iter + 1 << "]" << std::endl;
        
        // 1. State Estimation (MAP)
        // Compute MAP estimate using the current lambdas
        MAP_estimation_GMRF();
                
        // 2. Hyperparameter Maximization (LML & Gradient)
        prevLML = currentLML;
        
        std::pair<double, Eigen::Vector4d> result = calculate_LML_gradient();
        currentLML = result.first;
        gradientVec = result.second;
        
        // Control gradient explosion
        for (int i = 0; i < 4; ++i) {
            if (std::abs(gradientVec(i)) > 1e2) {
                gradientVec(i) = (gradientVec(i) > 0 ? 1e2 : -1e2);
            }
        }

        if (verbose) {
            std::cerr << "  -> LML: " << currentLML << std::endl;
            std::cerr << "  -> Lambdas: " << lambdaVec.transpose() << std::endl;
            std::cerr << "  -> Gradient: " << gradientVec.transpose() << std::endl;
        }

        // 3. Check Convergence
        if (iter > 0 && std::abs(currentLML - prevLML) < LML_threshold) {
            if (verbose) std::cerr << "[LML] Optimization converged. Delta LML: " << (currentLML - prevLML) << std::endl;
            break;
        }
        
        // 4. Gradient Ascent Update (Simplified Optimizer - but slow convergence)
        //lambdaVec += learningRate * gradientVec;

        // 4b. Update with Momentum & Adaptive Learning Rate (Step-Halving)
        if (iter > 0 && currentLML < prevLML) {
            // If LML decreased (it shouldn't in ascent), we reduce learning rate and velocity
            learningRate *= 0.5;
            velocity *= 0.5;
            if (verbose) std::cerr << "  (Step-back: Reducing learning rate to " << learningRate << ")" << std::endl;
        }

        velocity = momentum * velocity + learningRate * gradientVec;
        lambdaVec += velocity;

        // Apply boundary constraints (lambdas must be positive)
        for (int i = 0; i < 4; ++i) lambdaVec(i) = std::max(1e-6, lambdaVec(i));

        // Update class members for next MAP step
        lambdaObservations = lambdaVec(0);
        lambdaPrior_reg = lambdaVec(1);
        lambdaPrior_flux_conservation = lambdaVec(2);
        lambdaPrior_obstacles = lambdaVec(3);
    }
}


double CGMRF_map::getLambdaValue(FactorType type) const {
    switch (type) {
        case FactorType::Regularization: return lambdaPrior_reg;
        case FactorType::FluxConservation: return lambdaPrior_flux_conservation;
        case FactorType::Obstacle: return lambdaPrior_obstacles;
        case FactorType::Observation: return lambdaObservations;
        default: return 0.0;
    }
}

/*---------------------------------------------------------------
                MAXIMUM A POSTERIORI Estimation GMRF
  ---------------------------------------------------------------*/
void CGMRF_map::MAP_estimation_GMRF()
{
    try
    {
        /*
        * J (Jacobian) The J matrix contains the dr/dm for every factor in the graph
        *              J is size (nFactors x NumNodes)
        *
        * Lambda (weights) Is the Diagonal information matrix (contains the weights for each factor)
        *              Lambda is size (nFactors x nFactors).
        *
        * Y (vector of observations) contains the values of observations, 0 for prior factors
        *              y is size (nFactors x 1)
        *
        * R (Residuals) Since our system is deterministic, the residuals do not
        *              need to be re-evaluated on each iteration (we only perform 1 iteration).
        *              Therefore, R = -y, since we ALWAYS start from a all 0 map state.
        *
        * H (Hessian) = J' * Lambda * J
        *               H is size (NumNodes x NumNodes)
        *
        * G (gradient) = J' * Lambda * R
        *               g is size (NumNodes x 1)
        */
        
        if (estimateTiming) meanTimer.start();

        // 1. Get current number of factors (nPriorFactors is constant, but nObsFactors is dynamic)
        nFactors = nPriorFactors + nObsFactors;

        // 2. Setup the prior part of Jacobian (fixed)
        std::vector<Eigen::Triplet<double>> J_temp;
        J_temp.reserve(J.size() + nObsFactors);
        std::copy(J.begin(), J.end(), back_inserter(J_temp));

        // 3. Setup Lambda matrix and data (y) vector
        std::vector<Eigen::Triplet<double>> Lambda_temp;
        Lambda_temp.reserve(nFactors);
        Eigen::VectorXd y_temp;
        y_temp.resize(nFactors);
        y_temp.fill(0.0);
        
        // LAMBDA PRIOR FACTORS
        for (std::vector<std::pair<size_t, FactorType>>::iterator it = factor_types.begin(); it != factor_types.end(); ++it)
        {
            double lambda_value = getLambdaValue(it->second);
            Eigen::Triplet<double> lambda_entry(it->first, it->first, lambda_value);
            Lambda_temp.push_back(lambda_entry);
        }

        // 4. OBSERVATION FACTORS
        size_t count = nPriorFactors;       // start after the already introduced prior factors
        for (std::vector<TobservationGMRF>::iterator ito = activeObs.begin(); ito != activeObs.end(); ++ito)
        {
            // Each observation translates to 2 factors (Wx,Wy)
            // Wx range [1,N]
            Eigen::Triplet<double> J_entry(count, ito->cell_idx, 1);            
            J_temp.push_back(J_entry);
            y_temp[count] = ito->windX;
            //Eigen::Triplet<double> lambda_entry(count, count, ito->lambda);
            Eigen::Triplet<double> lambda_entry(count, count, lambdaObservations); // Use fixed lambda for observations (param)
            Lambda_temp.push_back(lambda_entry);
            count++;

            // Wy range [N+1,2N]            
            Eigen::Triplet<double> J_entry2(count, ito->cell_idx + N, 1);                
            J_temp.push_back(J_entry2);
            y_temp[count] = ito->windY;
            //Eigen::Triplet<double> lambda_entry2(count, count, ito->lambda);
            Eigen::Triplet<double> lambda_entry2(count, count, lambdaObservations); // Use fixed lambda for observations (param)
            Lambda_temp.push_back(lambda_entry2);
            count++;
        }

        
        // 5. Build Matrices (J, J', Lambda(A), H, G)
        Jsparse.resize(nFactors, 2 * N); // declares a column-major sparse matrix type of float
        Jsparse.setFromTriplets(J_temp.begin(), J_temp.end());
        if (false)
            std::cerr <<  "          [GMRF] Jsparse is (" << Jsparse.rows() << "," << Jsparse.cols() << ")" << std::endl;

        Eigen::SparseMatrix<double> JsparseT = Jsparse.transpose();
        if (false)
            std::cerr <<  "          [GMRF] JsparseT is (" << JsparseT.rows() << "," << JsparseT.cols() << ")" << std::endl;

        Eigen::SparseMatrix<double> Asparse(nFactors, nFactors); // declares a column-major sparse matrix type of float
        Asparse.setFromTriplets(Lambda_temp.begin(), Lambda_temp.end());
        if (false)
            std::cerr <<  "          [GMRF] Asparse is (" << Asparse.rows() << "," << Asparse.cols() << ")";

        Hsparse.resize(2 * N, 2 * N);
        Hsparse = JsparseT * Asparse * Jsparse; // size(2*N,2*N);
        if (false)
            std::cerr <<  "          [GMRF] Hsparse is (" << Hsparse.rows() << "," << Hsparse.cols() << ")" << std::endl;

        Eigen::VectorXd G = JsparseT * Asparse * y_temp;  // size(2*N,1);
        if (false)
            std::cerr <<  "          [GMRF] G is (" << G.rows() << "," << G.cols() << ")" << std::endl;
        
        
        //----------
        // 6. SOLVE H * m = G
        //----------
        // We use a Cholesky Factorization of Hessian --> chol( P * H * inv(P) )
        solver.compute(Hsparse);    // Computes the sparse Cholesky decomposition
        m_MAP_sol = solver.solve(G);
        if (false)
            std::cerr <<  "[GMRF] system solved with solution size (" << m_MAP_sol.rows() << "," << m_MAP_sol.cols() << ")" << std::endl;

        if (estimateTiming){
            meanTimer.stop(); // Stop Timer for Mean computation
            auto mean_time_ms = meanTimer.getMeanTimeMs();
            std::cerr << "[GMRF] Mean value estimated in " << mean_time_ms << " milliseconds" << std::endl;
        }
        
        // 7. Update GMRF values from current solution
        for (size_t j = 0; j < m_map.size(); j++)
        {
            m_map[j].mean = m_MAP_sol(j);   // Not iterative! no need to increment previous state
            m_map[j].std = 0.0;             // Not estimated yet. We need inv(H) diagonal for that.
        }

        // Calculate the final residual vector (r = J*m - y)
        residual = Jsparse * m_MAP_sol - y_temp;
        std::cerr << "[GMRF] MAP estimated with residual norm " << residual.norm() << std::endl;
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-updateMapEstimation_GMRF] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
    }
}


void CGMRF_map::computeUncertainty_GMRF()
{
    // The Hessian H size is (2*N,2*N)
    // The Covariance Matrix C is the inverse of the Hessian: C = H^{-1}.
    // The uncertainties (variances) are the diagonal elements of C.
    try
    {
        if (estimateTiming) stdTimer.start(); // Start Timer for Uncertainty computation        
        
        size_t matrix_size = Hsparse.rows();
        if (Hsparse.cols() != matrix_size || matrix_size == 0)
        {
            std::cerr << "[GMRF-computeUncertainty_GMRF] Error: Hsparse is not a square matrix or is empty." << std::endl;
            return;
        }

        if (verbose)
             std::cerr << "[GMRF] Computing uncertainty for matrix H of size (" << matrix_size << "," << matrix_size << ")..." << std::endl;

        // Use the same Cholesky decomposition used for the MAP estimation
        if (solver.info() != Eigen::Success)
        {
            // Failed to factorize the matrix. The system might be ill-conditioned or singular.
            std::cerr << "[GMRF-computeUncertainty_GMRF] Error: Failed to compute Cholesky factorization of Hsparse. Cannot compute uncertainty." << std::endl;

            if (estimateTiming)
                stdTimer.stop(); // Stop Timer for Uncertainty computation
            return;
        }

        // --- ITERATIVE METHOD ---
        // Not ideal, but... better than nothing
        // Vector for the right-hand side (RHS) of H * x_j = e_j
        Eigen::VectorXd e_j = Eigen::VectorXd::Zero(matrix_size);

        // 2. Loop through all state variables to solve H * x_j = e_j for each
        for (size_t j = 0; j < matrix_size; ++j)
        {
            // Set e_j to be the j-th standard basis vector (1 at j, 0 elsewhere)
            e_j.setZero();
            e_j(j) = 1.0;

            // Solve H * x_j = e_j. The j-th component of the solution x_j is the variance: (H^-1)_{j,j}.
            Eigen::VectorXd x_j = solver.solve(e_j);
            
            double variance = x_j(j);
            
            // 3. Update the standard deviation (uncertainty) in the map state
            if (variance > 0.0) {
                m_map[j].std = std::sqrt(variance); 
            } else {
                // Should not happen for a positive definite matrix, but handles numerical safety
                m_map[j].std = 0.0;
            }
        }
        
        if (verbose)
             std::cerr << "[GMRF] Uncertainty computation complete." << std::endl;

        if (estimateTiming){
            stdTimer.stop(); // Stop Timer for Uncertainty computation
            auto std_time_ms = stdTimer.getMeanTimeMs();
            std::cerr << "[GMRF] Uncertainty value estimated in " << std_time_ms << " milliseconds" << std::endl;
        }
    }
    catch (std::exception& e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-computeUncertainty_GMRF] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
    }
}


// ----------------------------------------------------------------------
// -------------------- LML & GRADIENT CALCULATION ----------------------
// ----------------------------------------------------------------------
std::pair<double, Eigen::Vector4d> CGMRF_map::calculate_LML_gradient()
{
    // LML = Term1 (Fit) + Term2 (Complexity) + Term3 (Prior Certainty)
    // ---------------------------------------------------------------
    double T1, T2, T3;
    Eigen::Vector4d gradientVec = Eigen::Vector4d::Zero();
    Eigen::Vector4d lambdaVec;
    lambdaVec << lambdaObservations, lambdaPrior_reg, lambdaPrior_flux_conservation, lambdaPrior_obstacles;

    // Compute T1 (Minimum Weighted Energy) and its gradient T1_grad
    // ---------------------------------------------------------------
    // T1 = -0.5 * residual' * Lambda * residual
    // T1_grad_i = -0.5 * sum_k (residual_k^2 * E_i(k,k))
    
    // We need the squared residual for each constraint, grouped by lambda type
    double r_sq_obs = 0.0;
    double r_sq_reg = 0.0;
    double r_sq_flux = 0.0;
    double r_sq_obst = 0.0;
    
    size_t count = 0;
    // Prior factors (first nPriorFactors rows)
    for (const auto& entry : factor_types) 
    {
        double r_sq = residual(entry.first) * residual(entry.first);
        switch (entry.second) {
            case FactorType::Regularization: r_sq_reg += r_sq; break;
            case FactorType::FluxConservation: r_sq_flux += r_sq; break;
            case FactorType::Obstacle: r_sq_obst += r_sq; break;
            default: break; // Should not happen
        }
        count = entry.first;
    }
    // Observation factors (rows nPriorFactors to nFactors-1)
    for (size_t i = nPriorFactors; i < nFactors; ++i) {
        r_sq_obs += residual(i) * residual(i);
    }
    
    // Calculate T1 value (Note: The Lambda multiplication is handled by the group sums)
    T1 = -0.5 * (lambdaObservations * r_sq_obs + lambdaPrior_reg * r_sq_reg + 
                 lambdaPrior_flux_conservation * r_sq_flux + lambdaPrior_obstacles * r_sq_obst);
                 
    // Calculate T1 gradient components (order: [0]Obs, [1]Reg, [2]Flux, [3]Obstacle)
    gradientVec(0) += -0.5 * r_sq_obs;
    gradientVec(1) += -0.5 * r_sq_reg;
    gradientVec(2) += -0.5 * r_sq_flux;
    gradientVec(3) += -0.5 * r_sq_obst;


    // Compute T2 (Complexity Penalty) and its gradient T2_grad
    // ---------------------------------------------------------------
    // T2 = -0.5 * log(|H|)
    // T2_grad_i = -0.5 * Tr(H^-1 * J' * E_i * J)

    // T2 Value: Use the Cholesky factorization (H = L L^T). |H| = prod(L_ii^2). log(|H|) = 2 * sum(log(L_ii)).
    double log_det_H_sum = 0.0;
    
    // Get the actual sparse matrix from the TriangularView (L factor)
    const auto& L_factor = solver.matrixL().nestedExpression();
    
    // Loop over the columns (or rows, since it's a view of a sparse matrix)
    for (int j = 0; j < L_factor.outerSize(); ++j) {
        // Iterate over the non-zero elements in the current column (j)
        for (Eigen::SparseMatrix<double>::InnerIterator it(L_factor, j); it; ++it) {
            if (it.row() == j) {
                // Found the diagonal element L_j,j. Log(|H|) = 2 * sum(log(L_j,j))
                log_det_H_sum += std::log(it.value());
                break; 
            }
        }
    }
    
    T2 = -1.0 * log_det_H_sum;  // Since log|H| = 2 * sum(log|Lii|)

    // T2 Gradient (Complexity Term) - Using Diagonal Approximation for Trace
    // Tr(H^-1 * J' * E_i * J) ~ sum_k ( (J' * E_i * J)_k,k * (H^-1)_k,k )
    // Since (H^-1)_k,k is H_diag_inv(k), we use that.
    
    // We need the diagonal of J' * E_i * J.
    // J' * E_i * J is a sparse matrix, so we compute the diagonal directly.
    
    Eigen::VectorXd J_E_J_diag = Eigen::VectorXd::Zero(Hsparse.rows());
    
    // The matrix J'*E_i*J is proportional to the Hessian H, but only using the factors weighted by lambda_i
    // J' * E_i * J = J' * (Lambda_i / lambda_i) * J = H_i / lambda_i
    
    // To implement J'*E_i*J for each i, we need to efficiently extract the Hessian contribution from J'*Lambda*J
    // The easiest way is to compute the diagonal of J'*E_i*J:
    // (J' * E_i * J)_k,k = sum_l (J_l,k^2 * E_i(l,l)) = sum_{l in constraints_i} J_l,k^2
    
    // Create a boolean mask for each factor type (e.g., E_obs, E_reg, etc.)
    std::vector<bool> is_obs_factor(nFactors, false);
    std::vector<bool> is_reg_factor(nFactors, false);
    std::vector<bool> is_flux_factor(nFactors, false);
    std::vector<bool> is_obst_factor(nFactors, false);

    for (const auto& entry : factor_types) {
        if (entry.second == FactorType::Regularization) is_reg_factor[entry.first] = true;
        if (entry.second == FactorType::FluxConservation) is_flux_factor[entry.first] = true;
        if (entry.second == FactorType::Obstacle) is_obst_factor[entry.first] = true;
    }
    for (size_t i = nPriorFactors; i < nFactors; ++i) {
        is_obs_factor[i] = true;
    }
    

    // Compute the diagonal of H^-1 (required for LML gradient approximation)
    Eigen::VectorXd H_diag_inv(Hsparse.rows());
    for (int j = 0; j < Hsparse.rows(); ++j) {
        Eigen::VectorXd e_j = Eigen::VectorXd::Zero(Hsparse.rows());
        e_j(j) = 1.0;
        H_diag_inv(j) = solver.solve(e_j)(j);
    }

    // Iterate over the columns of JsparseT to compute (J' * E_i * J)_k,k
    auto compute_trace_term = [&](const std::vector<bool>& factor_mask) -> double {
        double trace_term = 0.0;
        
        // Loop over the state dimensions (k)
        for (int k = 0; k < Hsparse.rows(); ++k) {
            // Compute (J' * E_i * J)_k,k = sum_{l in constraints_i} J_l,k^2
            double J_E_J_diag_k = 0.0;
            
            // Loop over J rows (l) for column k
            for (Eigen::SparseMatrix<double>::InnerIterator it(Jsparse, k); it; ++it) {
                size_t row_idx = it.row();
                if (row_idx < nFactors && factor_mask[row_idx]) {
                    J_E_J_diag_k += it.value() * it.value();
                }
            }
            
            // Apply the diagonal approximation: Tr(A H^-1) ~ sum_k A_k,k * (H^-1)_k,k
            trace_term += J_E_J_diag_k * H_diag_inv(k);
        }
        return trace_term;
    };
    
    double trace_obs  = compute_trace_term(is_obs_factor);
    double trace_reg  = compute_trace_term(is_reg_factor);
    double trace_flux = compute_trace_term(is_flux_factor);
    double trace_obst = compute_trace_term(is_obst_factor);
    
    // Apply Term 2 gradient: -0.5 * Tr(...)
    gradientVec(0) += -0.5 * trace_obs;
    gradientVec(1) += -0.5 * trace_reg;
    gradientVec(2) += -0.5 * trace_flux;
    gradientVec(3) += -0.5 * trace_obst;

    
    // Compute T3 (Prior Certainty) and its gradient T3_grad
    // ---------------------------------------------------------------
    // T3 = 0.5 * log(|Lambda|)
    // T3_grad_i = 0.5 * N_i / lambda_i
    
    // N_i: number of factors for each lambda type
    int N_obs = nObsFactors;
    int N_reg = 0;
    int N_flux = 0;
    int N_obst = 0;
    for (const auto& entry : factor_types) {
        if (entry.second == FactorType::Regularization) N_reg++;
        if (entry.second == FactorType::FluxConservation) N_flux++;
        if (entry.second == FactorType::Obstacle) N_obst++;
    }
    
    // T3 Value: log(|Lambda|) = sum_i (N_i * log(lambda_i))
    double log_det_Lambda = N_obs * std::log(lambdaObservations) + 
                            N_reg * std::log(lambdaPrior_reg) + 
                            N_flux * std::log(lambdaPrior_flux_conservation) + 
                            N_obst * std::log(lambdaPrior_obstacles);
    T3 = 0.5 * log_det_Lambda;
    
    // T3 Gradient
    gradientVec(0) += 0.5 * N_obs / lambdaObservations;
    gradientVec(1) += 0.5 * N_reg / lambdaPrior_reg;
    gradientVec(2) += 0.5 * N_flux / lambdaPrior_flux_conservation;
    gradientVec(3) += 0.5 * N_obst / lambdaPrior_obstacles;
    
    // Final LML (Term 1 + Term 2 + Term 3) - Ignore the constant C
    double LML = T1 + T2 + T3;

    return {LML, gradientVec};
}


WindVector CGMRF_map::getEstimation(int index)
{
    // GMRF stores the wind field as two separate variables per cell: Wx and Wy
    // index: cell index for Wx, index + N: cell index for Wy    
    double x = m_map[index].mean;
    double y = m_map[index + N].mean;
    double stdevX = std::max(0.001, m_map[index].std);
    double vx = stdevX * stdevX;
    double stdevY = std::max(0.001, m_map[index + N].std);
    double vy = stdevY * stdevY;

    // We convert to polar coordinates (module, direction) and propagate uncertainties
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan2(y, x);

    // Jacobian elements
    double dr_dx = x / r;
    double dr_dy = y / r;
    double dth_dx = -y / (r*r);
    double dth_dy = x / (r*r);

    // Uncertainty propagation: Sigma_polar = J * Sigma_cartesian * J^T
    // Assuming independence between x and y (no covariance) --> Sigma_cartesian = [vx 0; 0 vy]
    // Resulting Sigma_r^2 = (dr/dx)^2 * vx + (dr/dy)^2 * vy
    double var_r = (dr_dx * dr_dx * vx) + (dr_dy * dr_dy * vy);    
    // Resulting Sigma_theta^2 = (dth/dx)^2 * vx + (dth/dy)^2 * vy
    double var_theta = (dth_dx * dth_dx * vx) + (dth_dy * dth_dy * vy);

    return {
        r,                          // Module
        theta,                      // Direction
        sqrt(var_r),                // Sigma_mod
        sqrt(var_theta),            // Sigma_angle (radians)
        x,                          // Cartesian x
        y,                          // Cartesian y
        stdevX,                     // Sigma_x
        stdevY                      // Sigma_y
    };
}

WindVector CGMRF_map::getEstimation(double x, double y)
{
    int i = xy2idx(x, y);
    return getEstimation(i);
}


void CGMRF_map::save_grmf_factor_graph(std::vector<Eigen::Triplet<double>>& Jout, std::vector<Eigen::Triplet<double>>& Aout, Eigen::VectorXd& yout)
{
    bool save_dense = true;
    bool save_sparse = true;
    if (save_dense)
    {
        // 1. Jacobian
        Eigen::SparseMatrix<double> Jsparse(nFactors, 2 * N); // declares a column-major sparse matrix type of float
        Jsparse.setFromTriplets(Jout.begin(), Jout.end());
        Eigen::MatrixXd Jdense = Jsparse.toDense();

        // define the format you want, you only need one instance of this...
        const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
        // std::cerr <<  "[GMRF] Saving Factor-Graph to file...");
        std::ofstream file("~/gmrf_jacobian_dense.txt");
        if (file.is_open())
        {
            file << Jdense.format(CSVFormat) << '\n';
        }
        file.close();

        // 2. Information Matrix
        Eigen::SparseMatrix<double> Asparse(nFactors, nFactors); // declares a column-major sparse matrix type of float
        Asparse.setFromTriplets(Aout.begin(), Aout.end());
        Eigen::MatrixXd Adense = Asparse.toDense();
        std::ofstream file2("~/gmrf_lambda_dense.txt");
        if (file2.is_open())
        {
            file2 << Adense.format(CSVFormat) << '\n';
        }
        file2.close();
    }

    if (save_sparse)
    {
        std::cerr <<  "[GMRF] Saving Factor-Graph (list of triplets) to file..." << std::endl;
        std::cerr <<  "[GMRF] Jtriplets(" << Jout.size() << ",3), Atriplets(" << Aout.size() << ",3), numFactors(" << yout.rows() << ")" << std::endl;

        // 1. Jacobian
        std::ofstream file("~/gmrf_Jacobian.txt");
        if (file.is_open())
        {
            file << "# Jacobian of the GMRF: row col value"
                 << "\n";
            file << "# Num_cells_x = " << m_size_x << "\n";
            file << "# Num_cells_y = " << m_size_y << "\n";
            file << "# cell_size = " << m_resolution << "\n";
            for (std::vector<Eigen::Triplet<double>>::iterator it = Jout.begin(); it != Jout.end(); it++)
            {
                file << it->row() << " " << it->col() << " " << it->value() << "\n";
            }
        }
        file.close();

        // 2. Information Matrix
        std::ofstream file2("~/gmrf_Lambda.txt");
        if (file2.is_open())
        {
            for (std::vector<Eigen::Triplet<double>>::iterator it = Aout.begin(); it != Aout.end(); it++)
            {
                file2 << it->row() << " " << it->col() << " " << it->value() << "\n";
            }
        }
        file2.close();

        // 3. Save vector of observations
        std::ofstream file3("~/gmrf_observations.txt");
        if (file3.is_open())
        {
            file3 << yout;
        }
        file3.close();
    }
}

// Save Sparse matrices to file (for debug)
void CGMRF_map::save_grmf_factor_graph(Eigen::SparseMatrix<double>& H, Eigen::VectorXd& G)
{
    // 1. Hessian
    std::ofstream file("~/gmrf_hessian.txt");
    if (file.is_open())
    {
        file << "# Hessian of the GMRF: row col value"
             << "\n";
        for (int k = 0; k < H.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(H, k); it; ++it)
            {
                file << it.row() << " "; // row index
                file << it.col() << " "; // col index (here it is equal to k)
                file << it.value() << "\n";
            }
        }
    }
    file.close();

    // Hdense
    Eigen::MatrixXd Hdense = H.toDense();
    // define the format you want, you only need one instance of this...
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file3("~/gmrf_hessian_dense.txt");
    if (file3.is_open())
    {
        file3 << Hdense.format(CSVFormat) << '\n';
    }
    file3.close();

    // 2. Gradient
    std::ofstream file2("~/gmrf_gradient.txt");
    if (file2.is_open())
    {
        file2 << G;
    }
    file2.close();
}


int CGMRF_map::xy2idx(float x, float y) const
{
    int x_idx = static_cast<int>((x - m_x_min) / m_resolution);
    int y_idx = static_cast<int>((y - m_y_min) / m_resolution);
    return x_idx + y_idx * m_size_x;
}

// Public wrapper to expose cell center coordinates for visualization utilities
void CGMRF_map::id2xy_public(size_t id, double& x, double& y)
{
    id2xy(id, x, y);
}
