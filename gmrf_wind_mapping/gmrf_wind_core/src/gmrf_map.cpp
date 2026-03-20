#include "gmrf_wind_core/gmrf_map.h"
#include <cstdio>  // Necesario para fprintf
#include <iomanip> // Necesario para std::setprecision
#include <iostream>
#include <queue>

using namespace gmrfw;

/*---------------------------------------------------------------
                        Constructor
  ---------------------------------------------------------------*/
CGMRF_map::CGMRF_map(const TOccupancyMap& oc_map,
                     const Parameters& params,
                     bool verbose,
                     bool estimateTiming = false)
{
    // Set Params
    this->factor_select[0] = params.lambdaPrior_advection > 0;         // Advection
    this->factor_select[1] = params.lambdaPrior_mass_conservation > 0; // Mass Conservation
    this->factor_select[2] = params.lambdaPrior_diffusion > 0;         // Diffusion
    this->factor_select[3] = params.lambdaPrior_obstacles > 0;         // Obstacles
    this->verbose = verbose;
    this->estimateTiming = estimateTiming;
    update_lambdas(params.lambdaPrior_advection, params.lambdaPrior_mass_conservation, params.lambdaPrior_diffusion, params.lambdaPrior_obstacles);

    try
    {
        // Copy params to internal variables
        m_Ocgridmap = oc_map;     // Occupancy gridMap ( from ROS2 MapServer or other sources)
        m_resolution = params.cell_size; // Desired resolution to build the GMRF (m)

        //---------------------------------------
        // 0. Set GMRF Dimensions and Parameters
        //---------------------------------------
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
        N = m_size_x * m_size_y; // Number of cells in the GMRF

        std::cerr << "[GMRF_MAP] Generating GMRF for 2D WIND estimation..." << std::endl;

        //---------------------------------------
        // 1. Init the m_map container (2N cells)
        //---------------------------------------
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
            std::cerr << "--------------------------------" << std::endl;
        }

        //-------------------------------
        // 2. Memory Reservation (seepUp)
        //-------------------------------
        // Approximate number of static factors
        nPriorFactors = 8 * N;
        nObsFactors = 0;
        nFactors = nPriorFactors + nObsFactors;
        if (verbose)
        {
            std::cerr << "[CGMRF] Reserving Memory for Prior-factors" << std::endl;
            // Inform about the factors that will be set according to the selected priors
            std::string factor_names[4] = {"Advection", "Mass Conservation", "Diffusion", "Obstacles"};
            std::cerr << "[CGMRF] Prior factors selected: ";
            for (int i = 0; i < 4; i++)
            {
                if (factor_select[i])
                {
                    std::cerr << factor_names[i] << " ";
                }
            }
            std::cerr << std::endl;
        }

        // Reserve memory for the Jacobian and Information matrix
        J.clear();
        J.reserve(10 * nFactors); // Each factor accounts for the connectivity of multiple nodes -->multiple entries
        Lambda.clear();
        Lambda.reserve(nFactors); // Diagonal -> 1 entry for each factor
        factor_types.clear();

        //-------------------------
        // 3. Prior Factors
        //-------------------------
        size_t count = 0;
        int num_neighbors = 8; // Whether to consider 4 or 8 neighbors when setting factors.

        for (size_t j = 0; j < N; j++)
        {
            // Get cell_x and cell_y in the GMRF representation
            size_t jx, jy;
            id2cellxy(j, jx, jy);


            // --------------------------
            // ---   REGULARIZATION   ---
            // --------------------------
            //if (!is_cell_free(j))
            {
                // Force Wind estimation to be 0 (Damping) for numericall stability
                {
                    // Wx(j) = 0
                    J.push_back(Eigen::Triplet<double>(count, j, 1.0));
                    factor_types.push_back({count, FactorType::Regularization, j, j});
                    count++;

                    // Wy(j) = 0
                    J.push_back(Eigen::Triplet<double>(count, j + N, 1.0));
                    factor_types.push_back({count, FactorType::Regularization, j, j});
                    count++;
                }
            }

            // --------------------------------
            // ---   3.0 ADVECTION PRIOR    ---
            // --------------------------------
            if (factor_select[0])
            {
                // Wind in a cell should be similar to its neighbors in the direction of the wind.
                // This is a dynamic prior, so the weight (lambda) will depend on the current estimation
                // of the wind direction and will be updated at each iteration.
                struct RegNeighbor
                {
                    int dx;
                    int dy;
                    FactorType type;
                    bool is_diag;
                };
                std::vector<RegNeighbor> reg_nbs = {
                    {1, 0, FactorType::AdvectionHoriz, false}, // E
                    {0, 1, FactorType::AdvectionVert, false}   // N
                };

                if (num_neighbors == 8)
                {
                    // Consider also diagonals
                    reg_nbs.push_back({1, 1, FactorType::AdvectionDiag1, true});  // NE
                    reg_nbs.push_back({-1, 1, FactorType::AdvectionDiag2, true}); // NW
                }

                for (const auto& nb : reg_nbs)
                {
                    int nx = (int)jx + nb.dx;
                    int ny = (int)jy + nb.dy;
                    if (nx < 0 || nx >= m_size_x || ny < 0 || ny >= m_size_y)
                        continue;

                    // Neighbor idx
                    size_t nj = cellxy2id(nx, ny);

                    // Only if BOTH cells are free
                    if (is_cell_free(j) && is_cell_free(nj))
                    {
                        bool can_connect = true;
                        // If diagonal, check the "elbow" cells to avoid advection through corners
                        if (nb.is_diag)
                        {
                            if (!is_cell_free(cellxy2id(jx, ny)) || !is_cell_free(cellxy2id(nx, jy)))
                            {
                                can_connect = false;
                            }
                        }

                        if (can_connect)
                        {
                            // Link Wx
                            J.push_back(Eigen::Triplet<double>(count, j, 1.0));
                            J.push_back(Eigen::Triplet<double>(count, nj, -1.0));
                            factor_types.push_back({count, nb.type, j, nj});
                            count++;
                            // Link Wy
                            J.push_back(Eigen::Triplet<double>(count, j + N, 1.0));
                            J.push_back(Eigen::Triplet<double>(count, nj + N, -1.0));
                            factor_types.push_back({count, nb.type, j, nj});
                            count++;
                        }
                    }
                }
            }

            // --------------------------------
            // --- 3.1 MASS CONSERVATION    ---
            // --------------------------------
            // Eight Neighbor based Divergence free field constraint
            if (factor_select[1] && is_cell_free(j) && jx > 0 && jx < m_size_x - 1 && jy > 0 && jy < m_size_y - 1)
            {
                // As soon as any of its 8 clossest neighbour cells is free, set the factor
                bool set_divergence_free = false;

                // N,S,E,W
                if (is_cell_free(j - 1))
                {
                    Eigen::Triplet<double> J_entry(count, j - 1, -1);
                    J.push_back(J_entry);
                    set_divergence_free = true;
                }
                if (is_cell_free(j + 1)) // E
                {
                    Eigen::Triplet<double> J_entry(count, j + 1, 1);
                    J.push_back(J_entry);
                    set_divergence_free = true;
                }
                if (is_cell_free(j - m_size_x)) // S
                {
                    Eigen::Triplet<double> J_entry(count, j + N - m_size_x, -1);
                    J.push_back(J_entry);
                    set_divergence_free = true;
                }
                if (is_cell_free(j + m_size_x)) // N
                {
                    Eigen::Triplet<double> J_entry(count, j + N + m_size_x, 1);
                    J.push_back(J_entry);
                    set_divergence_free = true;
                }

                // Diagonals
                if (is_cell_free(j + m_size_x - 1)) // NW
                {
                    Eigen::Triplet<double> J_entry(count, j + m_size_x - 1, -0.5);
                    Eigen::Triplet<double> J_entry2(count, j + m_size_x - 1 + N, 0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set_divergence_free = true;
                }
                if (is_cell_free(j + m_size_x + 1)) // NE
                {
                    Eigen::Triplet<double> J_entry(count, j + m_size_x + 1, 0.5);
                    Eigen::Triplet<double> J_entry2(count, j + m_size_x + 1 + N, 0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set_divergence_free = true;
                }
                if (is_cell_free(j - m_size_x + 1)) // SE
                {
                    Eigen::Triplet<double> J_entry(count, j - m_size_x + 1, 0.5);
                    Eigen::Triplet<double> J_entry2(count, j - m_size_x + 1 + N, -0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set_divergence_free = true;
                }
                if (is_cell_free(j - m_size_x - 1)) // SW
                {
                    Eigen::Triplet<double> J_entry(count, j - m_size_x - 1, -0.5);
                    Eigen::Triplet<double> J_entry2(count, j - m_size_x - 1 + N, -0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set_divergence_free = true;
                }

                // If factor is to be set...
                if (set_divergence_free)
                {
                    factor_types.push_back({count, FactorType::MassConservation, j, j});
                    count++;
                    set_divergence_free = false;
                }
            }

            // --------------------------------
            // ---  3.2  VISCOUS DIFFUSION  ---
            // --------------------------------
            if (factor_select[2])
            {
                // Tikhonov (L2) regularization: smoothness prior. Force neighboring cells to have similar wind values
                struct RegNeighbor
                {
                    int dx;
                    int dy;
                    FactorType type;
                    bool is_diag;
                };
                std::vector<RegNeighbor> reg_nbs = {
                    {1, 0, FactorType::Diffusion, false}, // E
                    {0, 1, FactorType::Diffusion, false}  // N
                };
                /*
                if (num_neighbors == 8)
                {
                    // Consider also diagonals
                    reg_nbs.push_back({1, 1, FactorType::Diffusion, true});    // NE
                    reg_nbs.push_back({-1, 1, FactorType::Diffusion, true});   // NW
                }
                */

                for (const auto& nb : reg_nbs)
                {
                    int nx = (int)jx + nb.dx;
                    int ny = (int)jy + nb.dy;
                    if (nx < 0 || nx >= m_size_x || ny < 0 || ny >= m_size_y)
                        continue;

                    // neighbor idx
                    size_t nj = cellxy2id(nx, ny);

                    // Only regularize if BOTH cells are free
                    if (is_cell_free(j) && is_cell_free(nj))
                    {
                        bool can_connect = true;
                        // If diagonal, check the "elbow" cells to avoid smoothing through corners
                        if (nb.is_diag)
                        {
                            if (!is_cell_free(cellxy2id(jx, ny)) || !is_cell_free(cellxy2id(nx, jy)))
                            {
                                can_connect = false;
                            }
                        }

                        if (can_connect)
                        {
                            // Smooth Wx
                            J.push_back(Eigen::Triplet<double>(count, j, 1.0));
                            J.push_back(Eigen::Triplet<double>(count, nj, -1.0));
                            factor_types.push_back({count, nb.type, j, nj});
                            count++;
                            // Smooth Wy
                            J.push_back(Eigen::Triplet<double>(count, j + N, 1.0));
                            J.push_back(Eigen::Triplet<double>(count, nj + N, -1.0));
                            factor_types.push_back({count, nb.type, j, nj});
                            count++;
                        }
                    }
                }
            }

            // ----------------------
            // --- 3.3 OBSTACLES  ---
            // ----------------------
            if (factor_select[3])
            {
                // Goal: Force wind to zero at boundaries to prevent "bleeding" into walls.
                int dx_card[] = {1, 0}; // E, N
                int dy_card[] = {0, 1}; // E, N

                // Only consider E and N neighbors (never diagonals)
                for (int i = 0; i < 2; i++)
                {
                    // neighbor idx
                    int nx = (int)jx + dx_card[i];
                    int ny = (int)jy + dy_card[i];
                    if (nx < 0 || nx >= m_size_x || ny < 0 || ny >= m_size_y)
                        continue;

                    size_t nj = cellxy2id(nx, ny);

                    // If one is free and the neighbor is an obstacle (or vice-versa)
                    if (is_cell_free(j) != is_cell_free(nj))
                    {
                        size_t free_cell = is_cell_free(j) ? j : nj;

                        if (dx_card[i] != 0)
                        {
                            // Horizontal neighbor: Obstacle is to the side. Force Wx = 0
                            J.push_back(Eigen::Triplet<double>(count, free_cell, 1.0));
                            factor_types.push_back({count, FactorType::Obstacle, free_cell, free_cell});
                            count++;
                        }
                        else
                        {
                            // Vertical neighbor: Obstacle is above/below. Force Wy = 0
                            J.push_back(Eigen::Triplet<double>(count, free_cell + N, 1.0));
                            factor_types.push_back({count, FactorType::Obstacle, free_cell, free_cell});
                            count++;
                        }
                    }
                    else if (is_cell_free(j) && is_cell_free(nj) && !check_connectivity_between2cells(j, nj))
                    {
                        // An obstacle is in between cells j and j+1 (not seen due to resolution).
                        if (dx_card[i] != 0)
                        {
                            // Horizontal neighbor: Obstacle is to the side. Force Wx = 0
                            J.push_back(Eigen::Triplet<double>(count, j, 1.0));
                            factor_types.push_back({count, FactorType::Obstacle, j, j});
                            count++;
                            J.push_back(Eigen::Triplet<double>(count, nj, 1.0));
                            factor_types.push_back({count, FactorType::Obstacle, nj, nj});
                            count++;
                        }
                        else
                        {
                            // Vertical neighbor: Obstacle is above/below. Force Wy = 0
                            J.push_back(Eigen::Triplet<double>(count, j + N, 1.0));
                            factor_types.push_back({count, FactorType::Obstacle, j, j});
                            count++;
                            J.push_back(Eigen::Triplet<double>(count, nj + N, 1.0));
                            factor_types.push_back({count, FactorType::Obstacle, nj, nj});
                            count++;
                        }
                    }
                }
            }

        } // end for setting factors

        // Total Number of Factors
        //-------------------------
        nPriorFactors = count;
        nFactors = nPriorFactors + nObsFactors;
        activeObs.clear();
        std::cerr << "[CGMRF] Initialization Complete: " << nFactors << " factors for a map size of 2N=" << m_map.size() << " nodes" << std::endl;

        // Compute Distance to the closest obstacle for each cell (for later use in the Lambda factors)
        computeDistanceTransform();

        // DEBUG: Save to file
        //-------------------
        /*
            Eigen::VectorXd y_empty;
            y_empty.resize(nFactors);
            y_empty.fill(0.0);
            save_grmf_factor_graph(J,Lambda,y_empty);
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

void CGMRF_map::update_lambdas(double m_lambdaPrior_adv,
                               double m_lambdaPrior_mass,
                               double m_lambdaPrior_diff,
                               double m_lambdaPrior_obstacles)
{
    // Update internal precision parameters
    lambdaPrior_advection = m_lambdaPrior_adv;
    lambdaPrior_mass_conservation = m_lambdaPrior_mass;
    lambdaPrior_diffusion = m_lambdaPrior_diff;
    lambdaPrior_obstacles = m_lambdaPrior_obstacles;
}

void CGMRF_map::read_lambdas(double& m_lambdaPrior_adv,
                             double& m_lambdaPrior_mass,
                             double& m_lambdaPrior_diff,
                             double& m_lambdaPrior_obstacles)
{
    m_lambdaPrior_adv = lambdaPrior_advection;
    m_lambdaPrior_mass = lambdaPrior_mass_conservation;
    m_lambdaPrior_diff = lambdaPrior_diffusion;
    m_lambdaPrior_obstacles = lambdaPrior_obstacles;
}

/*---------------------------------------------------------------
                        Cell index transformations
  ---------------------------------------------------------------*/
// Get cell(x,y) in the GMRF representation from the general index in the array
void CGMRF_map::id2cellxy(size_t id, size_t& cell_x, size_t& cell_y) const
{
    cell_x = id % m_size_x;
    cell_y = (size_t)floor(id / m_size_x);
}

size_t CGMRF_map::cellxy2id(size_t cell_x, size_t cell_y) const
{
    return cell_x + cell_y * m_size_x;
}

// Get pose x,y (in meters) (in the GMRF representation) from the general index in the array
void CGMRF_map::id2xy(size_t id, double& x, double& y) const
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
bool CGMRF_map::is_cell_free(size_t id_gmrf) const
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

bool CGMRF_map::is_cell_boundary(size_t id_gmrf) const
{
    // Return true for cells at the borders of the map, which are considered as boundaries in CFD
    size_t cell_x, cell_y;
    id2cellxy(id_gmrf, cell_x, cell_y);
    return (cell_x == 0 || cell_x == m_size_x - 1 || cell_y == 0 || cell_y == m_size_y - 1);
}

/*---------------------------------------------------------------
             Check cell interconnectivity
  ---------------------------------------------------------------*/
// Check at OccupancyMap level, if two cells are interconnected, that is there exists no obstacles in between them.
bool CGMRF_map::check_connectivity_between2cells(size_t idx_1_gmrf, size_t idx_2_gmrf) const
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
        idx_1_oc = static_cast<int>((cell_1_x - m_Ocgridmap.origin_x) / m_Ocgridmap.resolution);                      // x component
        idx_1_oc += static_cast<int>((cell_1_y - m_Ocgridmap.origin_y) / m_Ocgridmap.resolution) * m_Ocgridmap.width; // y component
        idx_2_oc = static_cast<int>((cell_2_x - m_Ocgridmap.origin_x) / m_Ocgridmap.resolution);                      // x component
        idx_2_oc += static_cast<int>((cell_2_y - m_Ocgridmap.origin_y) / m_Ocgridmap.resolution) * m_Ocgridmap.width; // y component

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
bool CGMRF_map::insertObservation_GMRF(double wind_speed, double wind_direction, double var_wind_speed, double var_wind_direction, double x_pos, double y_pos)
{
    try
    {
        auto add_obs = [this](const TobservationGMRF& observation)
        {
            if (observation.cell_idx < 0 || observation.cell_idx > N)
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
            return false;

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
        new_obs.var_xx = var_r * cos(theta) * cos(theta) + r * r * var_theta * sin(theta) * sin(theta);
        new_obs.var_yy = var_r * sin(theta) * sin(theta) + r * r * var_theta * cos(theta) * cos(theta);
        new_obs.cov_xy = (var_r - r * r * var_theta) * sin(theta) * cos(theta);
        // Time variant or invariant observation
        new_obs.time_invariant = true; // Default behaviour, the obs will not lose weight with time.

        if (verbose)
            std::cerr << "[GMRF-MAP] New obs: Wx = " << new_obs.wind_x << " m/s Wy = " << new_obs.wind_y << " m/s" << std::endl;

        // Add Observation to GMRF
        add_obs(new_obs);
        nObsFactors += 2; // we add 2 factors foe each observation to account for Wx and Wy components
        return true;
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-insertObservation_GMRF] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
        return false;
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

void CGMRF_map::clearEstimation()
{
    try
    {
        // Clear any existing estimation data
        for (auto& cell : m_map)
        {
            cell.mean = 0.0;
            cell.var = 0.0;
            cell.covariance = 0.0;
        }
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[GMRF-clearEstimation] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
    }
}

void CGMRF_map::getObservationsIdx(std::vector<int>& obs_idx)
{
    obs_idx.clear();
    for (const auto& obs : activeObs)
    {
        obs_idx.push_back(obs.cell_idx);
    }
}


double CGMRF_map::getLambdaValue(FactorType type, size_t cell_idx, size_t neighbor_cell_idx, double ux, double uy, double neighbor_ux, double neighbor_uy) const
{
    // Advection weight for single cell
    auto calc_single_advective_weight = [&](size_t c_idx, double u, double v, double axis_x, double axis_y) 
    {
        // Wind magnitude
        double m_sq = u * u + v * v;

        // 1. Filter low wind magnitudes. No wind -> No advection
        if (std::isnan(u) || std::isnan(v) || m_sq < 1e-6) {
            return 0.0;
        }

        // 2. Magnitude and directional alignment (cosine of the angle between wind and axis)
        double m = std::sqrt(m_sq);
        double dot_product = (u * axis_x + v * axis_y);
        double alignment = std::abs(dot_product) / m; // [0,1], where 1 means perfectly aligned and 0 means perpendicular]

        // 0.9 = cos(25 degrees)
        if (alignment < 0.9) {
            return 0.0;
        }
        
        // 2. Smooth transition
        double dir_weight = std::pow(alignment, 8);
        
        // 3. RAYCASTING: Look for obstacles in the forward direction of the wind
        //----------------------------------------------------------
        int step_x = (axis_x > 0.1) ? 1 : ((axis_x < -0.1) ? -1 : 0);
        int step_y = (axis_y > 0.1) ? 1 : ((axis_y < -0.1) ? -1 : 0);

        if (dot_product < 0) {
            step_x = -step_x;
            step_y = -step_y;
        }

        size_t cx, cy;
        id2cellxy(c_idx, cx, cy);

        int dist_fwd = 0;
        int max_lookahead = 4;

        for (int s = 1; s <= max_lookahead; s++) 
        {
            int nx = (int)cx + s * step_x;
            int ny = (int)cy + s * step_y;

            // Check bounds
            if (nx < 0 || nx >= (int)m_size_x || ny < 0 || ny >= (int)m_size_y)
            {
                // Reached the edge of the map, but not obstacle (inlet/outlet)
                dist_fwd = max_lookahead; // No attenuation
                break;
            }

            if (!is_cell_free(cellxy2id(nx, ny))) break;
            
            dist_fwd++;
        }

        // 4. Spatial attenuation based on forward free space
        //if (dist_fwd < max_lookahead) 
        //    return 0.0; // No free space ahead -> no advection

        // Smooth transition: The more free space in the forward direction, the stronger the advection prior.
        double k = 0.8;
        double spatial_attenuation = 1.0 - std::exp(-k * dist_fwd);

        return dir_weight * spatial_attenuation;
    };
    
    auto calculate_symmetric_advection = [&](double axis_x, double axis_y) 
    {
        double weight_current = calc_single_advective_weight(cell_idx, ux, uy, axis_x, axis_y);
        double weight_neighbor = calc_single_advective_weight(neighbor_cell_idx, neighbor_ux, neighbor_uy, axis_x, axis_y);
        
        // Use max to account for upwind and downwind.
        return lambdaPrior_advection * std::max(weight_current, weight_neighbor);
    };

    switch (type)
    {
        case FactorType::AdvectionHoriz:
            return calculate_symmetric_advection(1.0, 0.0);
        case FactorType::AdvectionVert:
            return calculate_symmetric_advection(0.0, 1.0);
        case FactorType::AdvectionDiag1:
            return calculate_symmetric_advection(0.7071, 0.7071);
        case FactorType::AdvectionDiag2:
            return calculate_symmetric_advection(-0.7071, 0.7071);
            
        case FactorType::MassConservation: 
            return lambdaPrior_mass_conservation;
        case FactorType::Diffusion: 
            return lambdaPrior_diffusion;
        case FactorType::Obstacle: 
            return lambdaPrior_obstacles;
        case FactorType::Regularization:
            return 0.0001;
        default: 
        {
            std::cerr << "Warning: Unrecognized factor type in getLambdaValue." << std::endl;
            return 0.0001;
        }
    }
}


void CGMRF_map::computeDistanceTransform() 
{
    // This function computes the distance from each cell to the nearest obstacle using a breadth-first search (BFS) approach.

    m_cells_to_obs.assign(N, std::numeric_limits<int>::max());
    std::queue<size_t> q;

    // 1. Seed Obstacles (distance = 0)
    for (size_t j = 0; j < N; j++)
    {
        if (!is_cell_free(j))
        {
            m_cells_to_obs[j] = 0;
            q.push(j);
        }
    }

    // 2. 8-connectivity for Chebyshev steps
    int dx[] = {1, -1, 0, 0, 1, 1, -1, -1};
    int dy[] = {0, 0, 1, -1, 1, -1, 1, -1};

    while (!q.empty())
    {
        size_t curr_idx = q.front();
        q.pop();

        size_t cx, cy;
        id2cellxy(curr_idx, cx, cy);

        for (int i = 0; i < 8; i++)
        {
            int nx = (int)cx + dx[i];
            int ny = (int)cy + dy[i];

            // Check bounds
            if (nx >= 0 && nx < (int)m_size_x && ny >= 0 && ny < (int)m_size_y)
            {
                size_t next_idx = cellxy2id(nx, ny);
                if (m_cells_to_obs[next_idx] == std::numeric_limits<int>::max())
                {
                    m_cells_to_obs[next_idx] = m_cells_to_obs[curr_idx] + 1;
                    q.push(next_idx);
                }
            }
        }
    }
}

/*---------------------------------------------------------------
                MAXIMUM A POSTERIORI Estimation GMRF
  ---------------------------------------------------------------*/
void CGMRF_map::MAP_estimation_GMRF(int m_picard_iterations)
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
        if (verbose)
            std::cerr << "[CGMRF-MAP] Starting MAP Estimation with " << m_picard_iterations << " Picard-like iterations..." << std::endl;

        if (estimateTiming)
            meanTimer.start();

        // 1. Get current number of factors (nPriorFactors is constant, but nObsFactors is dynamic)
        nFactors = nPriorFactors + nObsFactors;

        // -----------------------------------------
        // ITERATIVE PICARD-LIKE APPROACH TO UPDATE LAMBDA
        // -----------------------------------------
        double prev_weighted_residual_norm = std::numeric_limits<double>::max();
        bool converged = false;
        for (int iter = 0; iter < m_picard_iterations; ++iter)
        {
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

            // LAMBDA PRIOR FACTORS (Dynamic values based on current wind estimation)
            for (std::vector<FactorInfo>::iterator it = factor_types.begin(); it != factor_types.end(); ++it)
            {
                // Get lambda value for the factor (type, cell_idx, and previous Wx,Wy wind estimation at cell_idx)          
                double lambda_value = getLambdaValue(it->type, it->cell_idx, it->neighbor_cell_idx, 
                    m_map[it->cell_idx].mean, m_map[it->cell_idx + N].mean,
                    m_map[it->neighbor_cell_idx].mean, m_map[it->neighbor_cell_idx + N].mean
                );
                Eigen::Triplet<double> lambda_entry(it->row_idx, it->row_idx, lambda_value);
                Lambda_temp.push_back(lambda_entry);
            }

            // 4. OBSERVATION FACTORS
            size_t count = nPriorFactors; // start after the already introduced prior factors
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
                    Eigen::Triplet<double> lambda_entry(count, count, 1.0 / ito->var_xx);
                    Lambda_temp.push_back(lambda_entry);
                    count++;

                    // Wy range [N+1,2N]
                    Eigen::Triplet<double> J_entry2(count, ito->cell_idx + N, 1);
                    J_temp.push_back(J_entry2);
                    y_temp[count] = ito->wind_y;
                    Eigen::Triplet<double> lambda_entry2(count, count, 1.0 / ito->var_yy);
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
                    if (det < 1e-9)
                        det = 1e-9; // Add small epsilon for numerical stability if r approx 0

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
                    Lambda_temp.push_back(Eigen::Triplet<double>(count, count, L_xx));         // diagonal x
                    Lambda_temp.push_back(Eigen::Triplet<double>(count + 1, count + 1, L_yy)); // diagonal y
                    Lambda_temp.push_back(Eigen::Triplet<double>(count, count + 1, L_xy));     // cross-term
                    Lambda_temp.push_back(Eigen::Triplet<double>(count + 1, count, L_xy));     // cross-term (symmetric)

                    count += 2;
                }
            }

            // 5. Build Matrices (J, J', Lambda(A), H, G)
            Jsparse.resize(nFactors, 2 * N); // declares a column-major sparse matrix type of float
            Jsparse.setFromTriplets(J_temp.begin(), J_temp.end());
            if (verbose)
                std::cerr << "          [GMRF] Jsparse is (" << Jsparse.rows() << "," << Jsparse.cols() << ")" << std::endl;

            Eigen::SparseMatrix<double> JsparseT = Jsparse.transpose();
            if (verbose)
                std::cerr << "          [GMRF] JsparseT is (" << JsparseT.rows() << "," << JsparseT.cols() << ")" << std::endl;

            Eigen::SparseMatrix<double> Asparse(nFactors, nFactors); // declares a column-major sparse matrix type of float
            Asparse.setFromTriplets(Lambda_temp.begin(), Lambda_temp.end());
            if (verbose)
                std::cerr << "          [GMRF] Asparse is (" << Asparse.rows() << "," << Asparse.cols() << ")" << std::endl;

            Hsparse.resize(2 * N, 2 * N);
            Hsparse = JsparseT * Asparse * Jsparse; // size(2*N,2*N);
            if (verbose)
                std::cerr << "          [GMRF] Hsparse is (" << Hsparse.rows() << "," << Hsparse.cols() << ")" << std::endl;

            Eigen::VectorXd G = JsparseT * Asparse * y_temp; // size(2*N,1);
            if (verbose)
                std::cerr << "          [GMRF] G is (" << G.rows() << "," << G.cols() << ")" << std::endl;

            //----------------------
            // 6. SOLVE H * m = G
            //----------------------
            // Sanity check for NaNs or Infs in Hsparse before factorization
            if (!Hsparse.coeffs().allFinite())
            {
                std::cerr << "=============================================================" << std::endl;
                std::cerr << "\033[1;31m" << "[GMRF-MAP] ERROR: Hessian matrix contains NaNs or Infs." << "\033[0m" << std::endl;
                std::cerr << "=============================================================" << std::endl;
                return;
            }

            // DEBUG: NUMERICAL STABILITY CHECK
            // ----------------------------------
            if (false)
            {
                // Check the diagonal values of Hsparse to identify potential zero or near-zero entries that could indicate singularity or ill-conditioning.
                std::cerr << "H max diag: " << Hsparse.diagonal().maxCoeff() << std::endl;
                std::cerr << "H min diag: " << Hsparse.diagonal().minCoeff() << std::endl;

                Eigen::VectorXd diag = Hsparse.diagonal();
                if (diag.size() != 2 * N)
                {
                    std::cerr << "Error de dimensiones: Diag " << diag.size() << " vs 2*N " << 2 * N << std::endl;
                    return;
                }

                int zero_count = 0;
                for (int i = 0; i < diag.size(); ++i)
                {
                    if (std::abs(diag(i)) < 1e-12)
                    {
                        zero_count++;
                        if (zero_count < 10)
                        { // Solo imprimimos los primeros 10 para no colapsar
                            int cell_idx = (i < N) ? i : i - N;
                            const char* comp = (i < N) ? "Wx" : "Wy";
                            std::cerr << "[GMRF] Nodo " << i << " (Celda " << cell_idx << " " << comp << ") es CERO." << std::endl;
                        }
                    }
                }
                std::cerr << "[GMRF] Total de nodos sin restricciones: " << zero_count << " de " << diag.size() << std::endl;

                // Add small value to diagonal for numerical stability (Tikhonov regularization)
                double epsilon = 1e-6;
                Hsparse.diagonal().array() += epsilon;
            }

            // Cholesky Factorization of Hessian (H = L * D * L')
            solver.compute(Hsparse); // Computes the sparse Cholesky decomposition
            if (solver.info() != Eigen::Success)
            {
                // In red color for better visibility
                std::cerr << "=============================================================" << std::endl;
                std::cerr << "\033[1;31m" << "[GMRF-MAP] ERROR: Failed to compute Cholesky factorization of Hsparse. The system might be ill-conditioned or singular. MAP estimation cannot proceed." << "\033[0m" << std::endl;
                std::cerr << "=============================================================" << std::endl;

                std::string ComputationInfo[4] = {"Success", "NumericalIssue", "NoConvergence", "InvalidInput"};
                std::cerr << "Eigen Computation Info Code: " << ComputationInfo[solver.info()] << std::endl;
                return;
            }

            // Get MAP solution: m_MAP_sol = H^-1 * G
            m_MAP_sol = solver.solve(G);
            if (false)
                std::cerr << "[GMRF] system solved with solution size (" << m_MAP_sol.rows() << "," << m_MAP_sol.cols() << ")" << std::endl;

            // 7. Update GMRF values from current iteration solution (MAP estimation)
            for (size_t j = 0; j < m_map.size(); j++)
            {
                // Under-relaxation update to improve convergence stability: new_mean = alpha * new_estimate + (1-alpha) * old_mean
                double alpha = 0.5; // Relaxation factor (0 < alpha <= 1)
                m_map[j].mean = alpha * m_MAP_sol(j) + (1 - alpha) * m_map[j].mean;
            }

            // 8. Residual Vector (r = J*m - y)
            residual = Jsparse * m_MAP_sol - y_temp;
            // Calculate the norm of the residual with the Asparse weights: sqrt(r' * A * r)
            double weighted_residual_norm = sqrt(residual.transpose() * Asparse * residual);

            if (verbose)
                std::cerr << "[GMRF] MAP estimated --> Picard iteration: " << iter + 1 << ". Weighted residual norm: " << weighted_residual_norm << std::endl;

            // Check convergence based on relative change in weighted residual norm
            double rel_change = std::abs(prev_weighted_residual_norm - weighted_residual_norm) / prev_weighted_residual_norm;
            if (rel_change < 1e-3) {
                std::cerr << "--> MAP Converged at iteration " << iter << std::endl;
                converged = true;
                break;
            }
            prev_weighted_residual_norm = weighted_residual_norm;

        } // end of Picard iterations

        if (!converged) {
            std::cerr << "--> MAP did not converge after " << m_picard_iterations << " iterations." << std::endl;
            std::cerr << "Final weighted residual norm: " << prev_weighted_residual_norm << std::endl;
        }

        if (estimateTiming)
        {
            meanTimer.stop(); // Stop Timer for Mean computation
            auto mean_time_ms = meanTimer.getMeanTimeMs();
            std::cerr << "[GMRF] MAP Mean value estimated in " << mean_time_ms << " milliseconds" << std::endl;
        }
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
        if (estimateTiming)
            stdTimer.start(); // Start Timer for Uncertainty computation

        size_t matrix_size = Hsparse.rows();
        size_t N = matrix_size / 2; // Number of cells

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

            m_map[j].var = std::max(0.0, col_x(j));            // Variance for Wx is the j-th diagonal element of H^-1
            m_map[j].covariance = std::max(0.0, col_x(j + N)); // Cross-term links Wx and Wy

            // --- Step 2: Solve for Wy column ---
            e_j.setZero();
            e_j(j + N) = 1.0;
            Eigen::VectorXd col_y = solver.solve(e_j);

            m_map[j + N].var = std::max(0.0, col_y(j + N)); // Variance for Wy is the (j+N)-th diagonal element of H^-1
            m_map[j + N].covariance = col_x(j + N);         // Cross-term is the same as before (symmetric)
        }

        if (verbose)
            std::cerr << "[GMRF] Uncertainty computation complete." << std::endl;

        if (estimateTiming)
        {
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
    double r = sqrt(pow(x, 2) + pow(y, 2)) + 1e-6; // Add small epsilon for numerical stability if r approx 0
    double theta = atan2(y, x);                    // Direction in radians, range [-pi, pi]

    // Jacobian elements: Cartesian to Polar transformation
    double dr_dx = x / r;
    double dr_dy = y / r;
    double dth_dx = -y / (r * r);
    double dth_dy = x / (r * r);

    // Uncertainty propagation: Sigma_polar = J * Sigma_cartesian * J^T
    double var_r = (dr_dx * dr_dx * vx) + (dr_dy * dr_dy * vy);
    double var_theta = (dth_dx * dth_dx * vx) + (dth_dy * dth_dy * vy);
    double cov_r_theta = (dr_dx * dth_dx * vx) + (dr_dy * dth_dy * vy) + (dr_dx * dth_dy + dr_dy * dth_dx) * cov_xy;

    return {
        r,           // Module
        theta,       // Direction
        var_r,       // Sigma_mod
        var_theta,   // Sigma_angle (radians)
        cov_r_theta, // Covariance between module and angle
        x,           // Cartesian x
        y,           // Cartesian y
        vx,          // Sigma_x
        vy,          // Sigma_y
        cov_xy       // Covariance between x and y
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
        std::cerr << "[GMRF] Saving Factor-Graph (list of triplets) to file..." << std::endl;
        std::cerr << "[GMRF] Jtriplets(" << Jout.size() << ",3), Atriplets(" << Aout.size() << ",3), numFactors(" << yout.rows() << ")" << std::endl;

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
void CGMRF_map::id2xy_public(size_t id, double& x, double& y) const
{
    id2xy(id, x, y);
}
