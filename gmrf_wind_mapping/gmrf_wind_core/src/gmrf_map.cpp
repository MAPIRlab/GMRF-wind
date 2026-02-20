#include "gmrf_wind_core/gmrf_map.h"
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
        init_cell.var = 0.0;
        init_cell.covariance = 0.0;
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
                        count++;

                        // Wy range [N+1,2N]                        
                        Eigen::Triplet<double> J_entry3(count, j + N, 1);
                        Eigen::Triplet<double> J_entry4(count, size_t(j + N + 1), -1);
                        J.push_back(J_entry3);
                        J.push_back(J_entry4);
                        factor_types.push_back({count, FactorType::Regularization});
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
                        count++;
                        
                        // Wx(j+1) = 0
                        Eigen::Triplet<double> J_entry2(count, j + 1, 1);
                        J.push_back(J_entry2);
                        factor_types.push_back({count, FactorType::Obstacle});
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
                }
                else if (is_cell_free(j + 1))
                {
                    // Cell j is occupied. Force Wx=0 at j+1
                    Eigen::Triplet<double> J_entry(count, j + 1, 1);
                    J.push_back(J_entry);
                    factor_types.push_back({count, FactorType::Obstacle});
                    count++;
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

                        // Wy range [N+1,2N]
                        Eigen::Triplet<double> J_entry3(count, j + N, 1);
                        Eigen::Triplet<double> J_entry4(count, j + N + m_size_x, -1);
                        J.push_back(J_entry3);
                        J.push_back(J_entry4);
                        factor_types.push_back({count, FactorType::Regularization});
                        count++;
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
                        
                        // Force Wy=0 at j+m_size_x
                        Eigen::Triplet<double> J_entry2(count, j + N + m_size_x, 1);
                        J.push_back(J_entry2);
                        factor_types.push_back({count, FactorType::Obstacle});
                        count++;
                    }
                }
                else if (is_cell_free(j))
                {
                    // Cell j+m_size_x is occupied. Force Wy=0 at j
                    Eigen::Triplet<double> J_entry(count, j + N, 1);                    
                    J.push_back(J_entry);
                    factor_types.push_back({count, FactorType::Obstacle});
                    count++;
                }
                else if (is_cell_free(j + m_size_x))
                {
                    // Cell j is occupied. Force Wy=0 at j+m_size_x
                    Eigen::Triplet<double> J_entry2(count, j + N + m_size_x, 1);
                    J.push_back(J_entry2);
                    factor_types.push_back({count, FactorType::Obstacle});
                    count++;
                }
                // else --> Both cells occupied -> Do nothing!
            }


            // 3.3 Factors for Law of Mass Conservation (avoid borders of map)
            //----------------------------------------------------------------------------------------
            if (is_cell_free(j) && jx > 0 && jx < m_size_x - 1 && jy > 0 && jy < m_size_y - 1)
            {
                // As soon as any of its 8 clossest neighbour cells is free, set the factor
                bool set = false;
                if (is_cell_free(j - 1)) // W
                {
                    Eigen::Triplet<double> J_entry(count, j - 1, -1);
                    J.push_back(J_entry);
                    set = true;
                }
                if (is_cell_free(j + 1)) // E
                {
                    Eigen::Triplet<double> J_entry(count, j + 1, 1);
                    J.push_back(J_entry);
                    set = true;
                }
                if (is_cell_free(j - m_size_x)) // S
                {
                    Eigen::Triplet<double> J_entry(count, j + N - m_size_x, -1);
                    J.push_back(J_entry);
                    set = true;
                }
                if (is_cell_free(j + m_size_x)) // N
                {
                    Eigen::Triplet<double> J_entry(count, j + N + m_size_x, 1);
                    J.push_back(J_entry);
                    set = true;
                }

                // Diagonals
                if (is_cell_free(j + m_size_x - 1)) // NW
                {
                    Eigen::Triplet<double> J_entry(count, j + m_size_x - 1, -0.5);
                    Eigen::Triplet<double> J_entry2(count, j + m_size_x - 1 + N, 0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }
                if (is_cell_free(j + m_size_x + 1)) // NE
                {
                    Eigen::Triplet<double> J_entry(count, j + m_size_x + 1, 0.5);
                    Eigen::Triplet<double> J_entry2(count, j + m_size_x + 1 + N, 0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }
                if (is_cell_free(j - m_size_x + 1)) // SE
                {
                    Eigen::Triplet<double> J_entry(count, j - m_size_x + 1, 0.5);
                    Eigen::Triplet<double> J_entry2(count, j - m_size_x + 1 + N, -0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }
                if (is_cell_free(j - m_size_x - 1)) // SW
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
                        double m_lambdaPrior_obstacles
                        )
{
    // Update internal variables
    lambdaPrior_reg = m_lambdaPrior_reg;
    lambdaPrior_flux_conservation = m_lambdaPrior_flux_conservation;
    lambdaPrior_obstacles = m_lambdaPrior_obstacles;
}

void CGMRF_map::read_lambdas(double &m_lambdaPrior_reg, double &m_lambdaPrior_flux_conservation, double &m_lambdaPrior_obstacles)
{
    m_lambdaPrior_reg = lambdaPrior_reg;
    m_lambdaPrior_flux_conservation = lambdaPrior_flux_conservation;
    m_lambdaPrior_obstacles = lambdaPrior_obstacles;
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


    try
    {
        // Check occupancy
        if (m_Ocgridmap.data[id_oc] >= 50.0)
            return false;
        else
            return true;
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
void CGMRF_map::insertObservation_GMRF(double wind_speed, double wind_direction, double var_wind_speed, double var_wind_direction, double x_pos, double y_pos)
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
        new_obs.wind_module = wind_speed;
        new_obs.wind_direction = wind_direction;
        new_obs.var_module = var_wind_speed;
        new_obs.var_direction = var_wind_direction;
        // Convert wind speed and direction to Wx and Wy components in the map reference system
        new_obs.wind_x = wind_speed * cos(wind_direction);
        new_obs.wind_y = wind_speed * sin(wind_direction);
        // Compute Full Covariance Sigma_xy (non-diagonal) from var_module and var_direction using error propagation
        double r = wind_speed;
        double theta = wind_direction;
        double var_r = var_wind_speed;
        double var_theta = var_wind_direction;
        new_obs.var_xx = var_r * cos(theta)*cos(theta) + r*r * var_theta * sin(theta)*sin(theta);
        new_obs.var_yy = var_r * sin(theta)*sin(theta) + r*r * var_theta * cos(theta)*cos(theta);
        new_obs.cov_xy = (var_r - r*r * var_theta) * sin(theta) * cos(theta);
        // Time variant or invariant observation
        new_obs.time_invariant = true; // Default behaviour, the obs will not lose weight with time.

        if (verbose)
            std::cerr << "[GMRF-MAP] New obs: Wx = " << new_obs.wind_x << " m/s Wy = " << new_obs.wind_y << " m/s" << std::endl;

        // Add Observation to GMRF
        add_obs(new_obs);
        nObsFactors += 2; // we add 2 factors foe each observation to account for Wx and Wy components
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-insertObservation_GMRF] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
    }
}

std::vector<TobservationGMRF> CGMRF_map::getObservations_GMRF()
{
    return activeObs;
}

void CGMRF_map::setObservations_GMRF(const std::vector<TobservationGMRF>& obs)
{
    activeObs.clear();
    activeObs.reserve(obs.size());
    for (const auto& element : obs)
        if (element.cell_idx >= 0 && element.cell_idx <= N)
            activeObs.push_back(element);
        else
            std::cerr << "[GMRF-MAP] Observation is outside of the map!" << std::endl;

    nObsFactors = activeObs.size() * 2;
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

double CGMRF_map::getLambdaValue(FactorType type) const {
    switch (type) {
        case FactorType::Regularization: return lambdaPrior_reg;
        case FactorType::FluxConservation: return lambdaPrior_flux_conservation;
        case FactorType::Obstacle: return lambdaPrior_obstacles;
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
            bool x_y_independent = false;

            if (x_y_independent)
            {
                // Each observation translates to 2 factors (Wx,Wy)
                // Wx range [1,N]
                Eigen::Triplet<double> J_entry(count, ito->cell_idx, 1);            
                J_temp.push_back(J_entry);
                y_temp[count] = ito->wind_x;
                Eigen::Triplet<double> lambda_entry(count, count, 1.0/ito->var_xx);
                Lambda_temp.push_back(lambda_entry);
                count++;

                // Wy range [N+1,2N]            
                Eigen::Triplet<double> J_entry2(count, ito->cell_idx + N, 1);                
                J_temp.push_back(J_entry2);
                y_temp[count] = ito->wind_y;
                Eigen::Triplet<double> lambda_entry2(count, count, 1.0/ito->var_yy);
                Lambda_temp.push_back(lambda_entry2);
                count++;
            }
            else
            {
                // If we consider correlation between Wx and Wy components of the same observation (derived from original Polar coordinates), 
                // we would need to introduce off-diagonal terms in Lambda.                    

                double var_xx = ito->var_xx;
                double var_yy = ito->var_yy;
                double cov_xy = ito->cov_xy;
                // 2. Invert 2x2 matrix to get Precision (Lambda)                
                double det = var_xx * var_yy - cov_xy * cov_xy;
                if (det < 1e-9) det = 1e-9; // Add small epsilon for numerical stability if r approx 0

                // Lambda_uv = inv(Sigma_xy) = (1/det) * [var_yy, -cov_xy; -cov_xy, var_xx]
                double L_xx = var_yy / det;
                double L_yy = var_xx / det;
                double L_xy = -cov_xy / det;

                // --- Update Triplets ---

                // Row for w_x=Obs_x component
                J_temp.push_back(Eigen::Triplet<double>(count, ito->cell_idx, 1));
                y_temp[count] = ito->wind_x;

                // Row for w_y=Obs_y component
                J_temp.push_back(Eigen::Triplet<double>(count + 1, ito->cell_idx + N, 1));
                y_temp[count + 1] = ito->wind_y;

                // Fill the 2x2 Precision Block in Lambda
                // This links row 'count' and 'count+1'
                Lambda_temp.push_back(Eigen::Triplet<double>(count,     count,     L_xx)); // diagonal x
                Lambda_temp.push_back(Eigen::Triplet<double>(count + 1, count + 1, L_yy)); // diagonal y
                Lambda_temp.push_back(Eigen::Triplet<double>(count,     count + 1, L_xy)); // cross-term
                Lambda_temp.push_back(Eigen::Triplet<double>(count + 1, count,     L_xy)); // cross-term (symmetric)

                count += 2;
            }
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
            m_map[j].var = 0.0;             // Not estimated yet.
            m_map[j].covariance = 0.0;             // Not estimated yet.
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
        size_t N = matrix_size / 2;     // Number of cells

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

        // 2x2 covariance blocks for Wx and Wy of each cell
        // ---------------------------------------------------------------------
        Eigen::VectorXd e_j = Eigen::VectorXd::Zero(matrix_size); // Standard basis vector for selecting columns of H^-1

        for (size_t j = 0; j < N; ++j)
        {
            // --- Step 1: Solve for Wx column ---
            e_j.setZero();
            e_j(j) = 1.0;
            Eigen::VectorXd col_x = solver.solve(e_j);
            
            m_map[j].var = std::max(0.0, col_x(j));             // Variance for Wx is the j-th diagonal element of H^-1
            m_map[j].covariance = std::max(0.0,col_x(j + N));   // Cross-term links Wx and Wy

            // --- Step 2: Solve for Wy column ---
            e_j.setZero();
            e_j(j + N) = 1.0;
            Eigen::VectorXd col_y = solver.solve(e_j);
            
            m_map[j+N].var = std::max(0.0, col_y(j + N));   // Variance for Wy is the (j+N)-th diagonal element of H^-1
            m_map[j+N].covariance = col_x(j + N);           // Cross-term is the same as before (symmetric)
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





WindVector CGMRF_map::getEstimation(int index)
{
    // GMRF stores the wind field as two separate variables per cell: Wx and Wy
    // index: cell index for Wx, index + N: cell index for Wy    
    double x = m_map[index].mean;
    double y = m_map[index + N].mean;
    double vx = std::max(0.0001, m_map[index].var);
    double vy = std::max(0.0001, m_map[index + N].var);
    double cov_xy = m_map[index].covariance; // covariance between Wx and Wy
    
    // We convert to polar coordinates (module, direction) and propagate uncertainties
    double r = sqrt(pow(x, 2) + pow(y, 2)) + 1e-6;  // Add small epsilon for numerical stability if r approx 0
    double theta = atan2(y, x);                     // Direction in radians, range [-pi, pi]

    // Jacobian elements: Cartesian to Polar transformation
    double dr_dx = x / r;
    double dr_dy = y / r;
    double dth_dx = -y / (r*r);
    double dth_dy = x / (r*r);

    // Uncertainty propagation: Sigma_polar = J * Sigma_cartesian * J^T
    double var_r = (dr_dx * dr_dx * vx) + (dr_dy * dr_dy * vy);
    double var_theta = (dth_dx * dth_dx * vx) + (dth_dy * dth_dy * vy);
    double cov_r_theta = (dr_dx * dth_dx * vx) + (dr_dy * dth_dy * vy) + (dr_dx * dth_dy + dr_dy * dth_dx) * cov_xy;

    return {
        r,                          // Module
        theta,                      // Direction
        var_r,                      // Sigma_mod
        var_theta,                  // Sigma_angle (radians)
        cov_r_theta,                // Covariance between module and angle
        x,                          // Cartesian x
        y,                          // Cartesian y
        vx,                         // Sigma_x
        vy,                         // Sigma_y
        cov_xy                      // Covariance between x and y
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
