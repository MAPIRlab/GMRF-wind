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
                     double m_lambdaPrior_flux_conservation, double m_lambdaPrior_obstacles,
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


/*---------------------------------------------------------------
                    updateMapEstimation_GMRF
  ---------------------------------------------------------------*/
void CGMRF_map::updateMapEstimation_GMRF(float lambdaObsLoss)
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

        // 2. Copy The prior part of Jacobian (fixed)
        std::vector<Eigen::Triplet<double>> J_temp;
        J_temp.reserve(J.size() + nObsFactors);
        std::copy(J.begin(), J.end(), back_inserter(J_temp));

        // 3. Generate the Lambda matrix with prior values of lambdas (parameters might have changed)
        std::vector<Eigen::Triplet<double>> Lambda_temp;
        Lambda_temp.reserve(factor_types.size() + nObsFactors);
        //std::copy(Lambda.begin(), Lambda.end(), back_inserter(Lambda_temp));
        for (std::vector<std::pair<size_t, FactorType>>::iterator it = factor_types.begin(); it != factor_types.end(); ++it)
        {
            double lambda_value = 0.0;
            switch (it->second)
            {
                case FactorType::Regularization:
                    lambda_value = lambdaPrior_reg;
                    break;
                case FactorType::FluxConservation:
                    lambda_value = lambdaPrior_flux_conservation;
                    break;
                case FactorType::Obstacle:
                    lambda_value = lambdaPrior_obstacles;
                    break;
                default:
                    lambda_value = 0.0;
                    break;
            }
            Eigen::Triplet<double> lambda_entry(it->first, it->first, lambda_value);
            Lambda_temp.push_back(lambda_entry);
        }

        // 4. Include Observations (Jacobian and Lambda)
        Eigen::VectorXd y_temp;
        y_temp.resize(nFactors);
        y_temp.fill(0.0);
        size_t count = nPriorFactors; // start after the already introduced prior factors
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

        // DEBUG - Save to file
        // save_grmf_factor_graph(J_temp,Lambda_temp,y_temp);

        // 5. Build Matrices (J, J', A, H, G)
        Eigen::SparseMatrix<double> Jsparse(nFactors, 2 * N); // declares a column-major sparse matrix type of float
        Jsparse.setFromTriplets(J_temp.begin(), J_temp.end());
        if (verbose)
            std::cerr <<  "          [GMRF] Jsparse is (" << Jsparse.rows() << "," << Jsparse.cols() << ")" << std::endl;

        Eigen::SparseMatrix<double> JsparseT; // size(2*N,nFactors);
        JsparseT = Eigen::SparseMatrix<double>(Jsparse.transpose());
        if (verbose)
            std::cerr <<  "          [GMRF] JsparseT is (" << JsparseT.rows() << "," << JsparseT.cols() << ")" << std::endl;

        Eigen::SparseMatrix<double> Asparse(nFactors, nFactors); // declares a column-major sparse matrix type of float
        Asparse.setFromTriplets(Lambda_temp.begin(), Lambda_temp.end());
        if (verbose)
            std::cerr <<  "          [GMRF] Asparse is (" << Asparse.rows() << "," << Asparse.cols() << ")";

        Eigen::SparseMatrix<double> Hsparse; // size(2*N,2*N);
        Hsparse = JsparseT * Asparse * Jsparse;
        //Hsparse = JsparseT * Jsparse;   // Since we fused Lambda into J
        if (verbose)
            std::cerr <<  "          [GMRF] Hsparse is (" << Hsparse.rows() << "," << Hsparse.cols() << ")" << std::endl;

        Eigen::VectorXd G = JsparseT * Asparse * y_temp;
        //Eigen::VectorXd G = JsparseT * y_temp; // Since we fused Lambda into J
        if (verbose)
            std::cerr <<  "          [GMRF] G is (" << G.rows() << "," << G.cols() << ")" << std::endl;
        // DEBUG - Save to file
        // save_grmf_factor_graph(Hsparse, G);
        

        
        //----------
        // 4. SOLVE
        //----------
        // We need to solve: H * m = G
        // We use a Cholesky Factorization of Hessian --> chol( P * H * inv(P) )
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(Hsparse);    // Computes the sparse Cholesky decomposition
        Eigen::VectorXd m_MAP_sol = solver.solve(G);

        if (estimateTiming){
            meanTimer.stop(); // Stop Timer for Mean computation
            auto mean_time_ms = meanTimer.getMeanTimeMs();
            std::cerr << "[GMRF] Mean value estimated in " << mean_time_ms << " milliseconds" << std::endl;
        }

        if (verbose)
            std::cerr <<  "[GMRF] system solved with solution size (" << m_MAP_sol.rows() << "," << m_MAP_sol.cols() << ")" << std::endl;
        
        // DEBUG - Save to file
        //std::ofstream file("~/gmrf_solution.txt");
        //if (file.is_open())
        //{
        //    file << m_inc;
        //}
        //file.close();


        // 5. Update GMRF values from current solution
        for (size_t j = 0; j < m_map.size(); j++)
        {
            m_map[j].mean = m_MAP_sol(j);   // Not iterative! no need to increment previous state
            m_map[j].std = 0.0;             // Not estimated yet. We need inv(H) diagonal for that.
        }

        // Uncertainty estimation
        if (estimateTiming) stdTimer.start(); // Start Timer for Uncertainty computation        
        computeUncertainty_GMRF(Hsparse);
        if (estimateTiming){
            stdTimer.stop(); // Stop Timer for Uncertainty computation
            auto std_time_ms = stdTimer.getMeanTimeMs();
            std::cerr << "[GMRF] Uncertainty value estimated in " << std_time_ms << " milliseconds" << std::endl;
        }

        // 6. Update Information/Strength of Active Observations
        //-------------------------------------------------------
        std::vector<TobservationGMRF>::iterator ito = activeObs.begin();
        while (ito != activeObs.end())
        {
            if (ito->time_invariant == false)
            {
                ito->lambda -= lambdaObsLoss;
                if (ito->lambda <= 0.0)
                {
                    ito = activeObs.erase(ito);
                    nObsFactors -= 2;
                }
                else
                    ++ito;
            }
            else
                ++ito;
        }
        if (verbose)
            std::cerr <<  "[GMRF] "<< nObsFactors << " ObservationFactors are active" << std::endl;
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-updateMapEstimation_GMRF] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
    }
}


void CGMRF_map::computeUncertainty_GMRF(const Eigen::SparseMatrix<double>& Hsparse)
{
    // The Hessian H is the Information Matrix (Lambda) for the MAP estimate. Here, H size is (2*N,2*N)
    // The Covariance Matrix C is the inverse of the Hessian: C = H^{-1}.
    // The uncertainties (variances) are the diagonal elements of C.

    try
    {
        size_t matrix_size = Hsparse.rows();
        if (Hsparse.cols() != matrix_size || matrix_size == 0)
        {
            std::cerr << "[GMRF-computeUncertainty_GMRF] Error: Hsparse is not a square matrix or is empty." << std::endl;
            return;
        }

        if (verbose)
             std::cerr << "[GMRF] Computing uncertainty for matrix H of size (" << matrix_size << "," << matrix_size << ")..." << std::endl;

        // Use the same SimplicialLLT solver used for the mean computation.
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(Hsparse);

        if (solver.info() != Eigen::Success)
        {
            // Failed to factorize the matrix. The system might be ill-conditioned or singular.
            std::cerr << "[GMRF-computeUncertainty_GMRF] Error: Failed to compute Cholesky factorization of Hsparse. Cannot compute uncertainty." << std::endl;
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
    }
    catch (std::exception& e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-computeUncertainty_GMRF] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
    }
}


WindVector CGMRF_map::getEstimation(int index)
{
    double module = sqrt(pow(m_map[index].mean, 2) + pow(m_map[index + N].mean, 2));
    double direction = atan2(m_map[index + N].mean, m_map[index].mean);
    double stdev = std::max(0.1, sqrt(pow(m_map[index].std, 2) + pow(m_map[index + N].std, 2)));

    return {module, direction, stdev};
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
