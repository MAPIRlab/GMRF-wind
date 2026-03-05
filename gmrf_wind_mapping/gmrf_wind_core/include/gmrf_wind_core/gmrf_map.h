/*  ================================================================
    Core implementation of a GMRF (map)
    GMRF class implements the probability map and methods 
    for adapting to an occupancy gridmap, insterting new
    observations and update the map. It also includes methods
    for visualization and reading the map estimations at specified
    locations.
    ================================================================
*/
#pragma once
#include <eigen3/Eigen/Sparse>
#include <fstream> // std::ofstream
#include <math.h>  /* atan2 */
#include <numeric> // std::accumulate
#include <chrono>   // Needed for time measurements

// Data structure for each cell in the GMRF
struct TRandomFieldCell
{
    double mean; 
    double var;         // Variance of the estimation (std^2)
    double covariance;  // Covariance between Wx and Wy if needed
};

// Data structure to hold the occupancy map to which the GMRF is adapted
struct TOccupancyMap
{
    std::vector<int8_t> data;
    double resolution;      // cell size (m)
    size_t width;           // number of cells in X direction
    size_t height;          // number of cells in Y direction
    // Real-world pose of the cell (0,0) in the map ref system
    double origin_x;        // world coordinates of the origin of the map
    double origin_y;        // world coordinates of the origin of the map
};

// Data structure to return wind vector estimations
struct WindVector
{
    // polar coordinates
    double module;
    double angle; // in radians, in the map reference system (0 rad means wind blowing towards +X, pi/2 means wind blowing towards +Y)
    double varMod;
    double varAngle;
    double covModAngle;
    // Cartesian coordinates
    double x;
    double y;
    double varX;
    double varY;
    double covXY; // covariance between x and y components

    Eigen::Vector2d asEigen()
    {
        return Eigen::Vector2d(x, y);
    }
};


// Struct to hold timing data
struct TimeStats {
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::vector<long long> timings_ms;

    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    void stop() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        timings_ms.push_back(duration.count());
    }

    double getMeanTimeMs() const {
        if (timings_ms.empty()) {
            return 0.0;
        }
        // Accumulate
        long long sum_ms = std::accumulate(timings_ms.begin(), timings_ms.end(), 0LL);
        return static_cast<double>(sum_ms) / timings_ms.size();
    }
};

// Enum for the Factor Types
enum class FactorType {
    Advection,
    AdvectionHoriz,
    AdvectionVert,
    AdvectionDiag1,    // +45 deg
    AdvectionDiag2,    // 135 deg
    MassConservation,
    Diffusion,
    Vorticity,
    Obstacle,
    Observation,
    //Regularization,
    //RegularizationHoriz,
    //RegularizationVert,
    //RegularizationDiag1,    // +45 deg
    //RegularizationDiag2     // 135 deg
};

struct FactorInfo {
    size_t row_idx;    // The 'count' to identify the row in J and Lambda
    FactorType type;   // Type of factor (MassConservation, Vorticity, etc.)
    size_t cell_idx;   // The 'j' index of the central cell associated with the factor
};

class CGMRF_map
{
public:
    // Create GMRF from an occupancy gridmap
    CGMRF_map(const TOccupancyMap& oc_map, 
                    float cell_size,
                    bool factor_select[5],
                    bool verbose,
                    bool estimateTiming);
    ~CGMRF_map();

    // Observations and Parameters
    void insertObservation_GMRF(double wind_speed, double wind_direction, double var_wind_speed, double var_wind_direction, double x_pos, double y_pos);
    void clearObservations_GMRF();
    void getObservationsIdx(std::vector<int>& obs_idx);
    void update_lambdas(double m_lambdaPrior_adv, double m_lambdaPrior_mass, double m_lambdaPrior_diff, double m_lambdaPrior_vorticity, double m_lambdaPrior_obstacles);
    void read_lambdas(double &m_lambdaPrior_adv, double &m_lambdaPrior_mass, double &m_lambdaPrior_diff, double &m_lambdaPrior_vorticity, double &m_lambdaPrior_obstacles);
   
    // GMRF Estimation
    void MAP_estimation_GMRF();         // Solves the Least Squares linear system (MAP estimator)
    void computeUncertainty_GMRF();

    // Read estimation
    WindVector getEstimation(int index);
    WindVector getEstimation(double x, double y);

    // Public accessor to get cell center coordinates (x,y) in meters from cell index
    // This forwards to the internal id2xy() utility and is provided for visualization
    // helpers that need world coordinates for each GMRF cell.
    void id2xy_public(size_t id, double& x, double& y);
    bool is_cell_free(size_t id_gmrf);
     
    Eigen::Vector2i map_size()
    {
        return {m_size_x, m_size_y};
    }

    float map_resolution()
    {
        return m_resolution;
    }

    std::array<double,4> map_dimensions_meters()
    {
        return {m_x_min, m_x_max, m_y_min, m_y_max};
    }

protected:
    std::vector<TRandomFieldCell> m_map;      // GMRF container of nodes
    TOccupancyMap m_Ocgridmap;                // Occupancy gridmap of the environment
    float m_x_min, m_x_max, m_y_min, m_y_max; // dimensions (m)
    float m_resolution;                       // cell_size (m)
    size_t m_size_x, m_size_y;                // dimensions in CellNumber
    size_t N;                                 // number of cells in the GMRF (we have 2N nodes)
    
    // Time Stats
    bool estimateTiming;
    TimeStats meanTimer;
    TimeStats stdTimer;

    // GMRF
    size_t nPriorFactors;                 // Static factors (dont change over time)
    size_t nObsFactors;                   // Dynamic factors due to new observations
    size_t nFactors;                      // Total num of factors
    bool verbose;
    bool factor_select[5];                // Whether to set or not the factors of each type (mass conservation, vorticity, obstacles, regularization, advection)
    double lambdaPrior_advection;         // Weight for advection prior (wind in a cell should be similar to its neighbors in the direction of the wind, this is dynamic as it depends on the current estimation of the wind direction)
    double lambdaPrior_mass_conservation; // Weight for mass conservation (divergence free)
    double lambdaPrior_diffusion;        // Weight for difussion prior (wind in a cell should be similar to its neighbors in all directions, also regularization)
    double lambdaPrior_vorticity;         // Weight for vorticity prior (curl free) This one is dynamic, so this is the max value.
    double lambdaPrior_obstacles;         // Weight for wind close to obstacles prior
    int m_picard_iterations = 3;          // Threshold/count for Picard Method
    std::vector<int> m_cells_to_obs;      // Distance in cells to the closest obstacle
    
    struct TobservationGMRF
    {
        size_t cell_idx;
        // Sensor data in polar coordinates (wind speed and direction)
        double wind_module;
        double wind_direction;
        double var_module;
        double var_direction;
        // Cartesian coordinates for the latent space (windX and windY) + full covariance
        double wind_x;
        double wind_y;
        double var_xx;
        double var_yy;
        double cov_xy;
        //double lambda;
        bool time_invariant; // if the observation will lose weight (lambda) as time goes on (default false)
    };

    // GMRF Matrices and Structures
    std::vector<Eigen::Triplet<double>> J;          // Jacobian
    std::vector<Eigen::Triplet<double>> Lambda;     // the information matrix (weights or inverse of covariances)
    std::vector<FactorInfo> factor_types;           // Vector with the type of each factor
    std::vector<TobservationGMRF> activeObs;        // Vector with the active observations and their respective Information
    Eigen::SparseMatrix<double> Jsparse;            // Sparse Jacobian matrix
    Eigen::SparseMatrix<double> Hsparse;            // Sparse Information matrix (H = J' * Lambda * J)
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;        // Cholesky solver for the linear system
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;       // LDLT Cholesky solver for the linear system (allows for non-positive definite matrices)
    Eigen::VectorXd gradient;                       // Gradient vector (G = J' * Lambda * y)
    Eigen::VectorXd m_MAP_sol;                     // Solution vector for the linear system (m_MAP)
    Eigen::VectorXd residual;                      // Residual vector (r = J*m_MAP_sol - y)

    // Util Functions
    bool check_connectivity_between2cells(size_t idx_1_gmrf, size_t idx_2_gmrf);
    double getLambdaValue(FactorType type, size_t cell_idx, double ux, double uy) const;
    void computeDistanceTransform();
    
    void id2cellxy(size_t id, size_t& cell_x, size_t& cell_y);
    size_t cellxy2id(size_t cell_x, size_t cell_y);
    int xy2idx(float x, float y) const;
    void id2xy(size_t id, double& x, double& y);

    // Visualization
    void save_grmf_factor_graph(std::vector<Eigen::Triplet<double>>& Jout, std::vector<Eigen::Triplet<double>>& Aout, Eigen::VectorXd& yout);
    void save_grmf_factor_graph(Eigen::SparseMatrix<double>& H, Eigen::VectorXd& G);
};