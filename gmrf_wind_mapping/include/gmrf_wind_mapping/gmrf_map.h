/*  ================================================================
    Core implementation of a GMRF (map)
    GMRF class implements the probability map and methods 
    for adapting to an occupancy gridmap, insterting new
    observations and update the map. It also includes methods
    for visualization and reading the map estimations at specified
    locations.
    ================================================================
*/
#include <eigen3/Eigen/Sparse>
#include <fstream> // std::ofstream
#include <math.h>  /* atan2 */
#include <numeric> // std::accumulate
#include <chrono>   // Needed for time measurements

// Data structure for each cell in the GMRF
// Stores mean and standard deviation, as we build a "Gaussian" Random Field
struct TRandomFieldCell
{
    double mean;
    double std;
};

// Data structure to hold the occupancy map to which the GMRF is adapted
struct TOccupancyMap
{
    std::vector<int8_t> data;
    double resolution;      // cell size (m)
    size_t width;        // number of cells in X direction
    size_t height;       // number of cells in Y direction
    // Real-world pose of the cell (0,0) in the map ref system
    double origin_x;        // world coordinates of the origin of the map
    double origin_y;        // world coordinates of the origin of the map
};

// Data structure to return wind vector estimations
struct WindVector
{
    // polar coordinates
    double module;
    double direction;
    double stdMod;
    double stdAngle;
    // Cartesian coordinates
    double x;
    double y;
    double stdX;
    double stdY;

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
    Regularization,
    Obstacle,
    FluxConservation,
    Observation
};

class CGMRF_map
{
public:
    // Create GMRF from an occupancy gridmap
    // And sets the prior weights for the different factors
    CGMRF_map(  const TOccupancyMap& oc_map, 
                float cell_size,
                double m_lambdaPrior_reg,
                double m_lambdaPrior_flux_conservation, 
                double m_lambdaPrior_obstacles,
                double m_lambdaObservations,
                bool verbose,
                bool estimateTiming
            );
    ~CGMRF_map();

    // Observations and Parameters
    void insertObservation_GMRF(double wind_speed, double wind_direction, double x_pos, double y_pos, double lambdaObs);
    void clearObservations_GMRF();
    void update_lambdas(double m_lambdaPrior_reg, double m_lambdaPrior_flux_conservation, double m_lambdaPrior_obstacles, double m_lambdaObservations);
    void read_lambdas(double &m_lambdaPrior_reg, double &m_lambdaPrior_flux_conservation, double &m_lambdaPrior_obstacles, double &m_lambdaObservations);
   
    // GMRF Estimation
    void optimize_LML(bool performLMLOptimization, int maxIterations, double learningRate, double LML_threshold); 
    void MAP_estimation_GMRF();         // Solves the Least Squares linear system (MAP estimator)
    void computeUncertainty_GMRF();
    std::pair<double, Eigen::Vector4d> calculate_LML_gradient();

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
    bool verbose;

    // Time Stats
    bool estimateTiming;
    TimeStats meanTimer;
    TimeStats stdTimer;

    // GMRF
    size_t nPriorFactors;                 // Static factors (dont change over time)
    size_t nObsFactors;                   // Dynamic factors due to new observations
    size_t nFactors;                      // Total num of factors
    double lambdaPrior_reg;               // Weight for regularization prior -> neighbour cells have similar wind vectors
    double lambdaPrior_flux_conservation; // Weight for flux conservation law prior
    double lambdaPrior_obstacles;         // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind
    double lambdaObservations;            // Weight for observations
    // SQRT values (to build J matrix without the Lambda diagonal matrix)
    //double lambdaPrior_reg_sqrt;               // Weight for regularization prior -> neighbour cells have similar wind vectors
    //double lambdaPrior_mass_conservation_sqrt; // Weight for mass conservation law prior
    //double lambdaPrior_obstacles_sqrt;         // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind

    struct TobservationGMRF
    {
        size_t cell_idx;
        double windX;
        double windY;
        double lambda;
        bool time_invariant; // if the observation will lose weight (lambda) as time goes on (default false)
    };

    // GMRF Matrices and Structures
    std::vector<Eigen::Triplet<double>> J;          // Jacobian
    std::vector<Eigen::Triplet<double>> Lambda;     // the information matrix (weights or inverse of covariances)
    std::vector<std::pair<size_t, FactorType>> factor_types; // Is the Lambda matrix, but with the type of each factor
    std::vector<TobservationGMRF> activeObs;        // Vector with the active observations and their respective Information
    Eigen::SparseMatrix<double> Jsparse;            // Sparse Jacobian matrix
    Eigen::SparseMatrix<double> Hsparse;            // Sparse Information matrix (H = J' * Lambda * J)
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;        // Cholesky solver for the linear system
    Eigen::VectorXd gradient;                       // Gradient vector (G = J' * Lambda * y)
    Eigen::VectorXd m_MAP_sol;                     // Solution vector for the linear system (m_MAP)
    Eigen::VectorXd residual;                      // Residual vector (r = J*m_MAP_sol - y)

    // Util Functions
    bool check_connectivity_between2cells(size_t idx_1_gmrf, size_t idx_2_gmrf);
    double getLambdaValue(FactorType type) const;

    int xy2idx(float x, float y) const;
    void id2cellxy(size_t id, size_t& cell_x, size_t& cell_y);
    void id2xy(size_t id, double& x, double& y);

    // Visualization
    void save_grmf_factor_graph(std::vector<Eigen::Triplet<double>>& Jout, std::vector<Eigen::Triplet<double>>& Aout, Eigen::VectorXd& yout);
    void save_grmf_factor_graph(Eigen::SparseMatrix<double>& H, Eigen::VectorXd& G);
};