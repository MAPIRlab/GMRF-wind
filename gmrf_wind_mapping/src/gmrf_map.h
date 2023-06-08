
#include "rclcpp/rclcpp.hpp"
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <eigen3/Eigen/Sparse>
#include <fstream>                              // std::ofstream
#include "visualization_msgs/msg/marker.hpp"
#include "visualization_msgs/msg/marker_array.hpp"
#include <math.h>                               /* atan2 */


#define NUM_CELL_TEMPLATES 200      //For plotting only
struct TRandomFieldCell
{
    double mean;
    double std;
};



/** GMRF class implementing the probability map and methods for insterting new observations and update the map */
class CGMRF_map
{
public:
    CGMRF_map(rclcpp::Node* _node, const nav_msgs::msg::OccupancyGrid &oc_map,
                         float cell_size,
                         double m_lambdaPrior_reg,
                         double m_lambdaPrior_mass_conservation,
                         double m_lambdaPrior_obstacles,
                         std::string m_colormap, int max_points_cell);
    ~CGMRF_map();

    //insert new observation
    void  insertObservation_GMRF(double wind_speed, double wind_direction, double x_pos, double y_pos, double lambdaObs);

    //solves the minimum quadratic system to determine the new concentration of each cell
    void  updateMapEstimation_GMRF(float lambdaObsLoss);

    //Visualization
    void get_as_markerArray(visualization_msgs::msg::MarkerArray &ma, std::string frame_id);

    Eigen::Vector3d getEstimation(double x, double y);
protected:
    rclcpp::Node* node;
    std::vector<TRandomFieldCell>           m_map;                                  // GMRF container of nodes
    nav_msgs::msg::OccupancyGrid                 m_Ocgridmap;                            // Occupancy gridmap of the environment
    float                                   m_x_min,m_x_max,m_y_min,m_y_max;        // dimensions (m)
    float                                   m_resolution;                           // cell_size (m)
    size_t                                  m_size_x, m_size_y;                     // dimensions in CellNumber
    size_t                                  N;                                      // number of cells in the GMRF (we have 2N nodes)
    float                                   res_coef;

    //GMRF
    size_t                                  nPriorFactors;                      // Static/fixed factors
    size_t                                  nObsFactors;                        // Dynamic factors due to observations
    size_t                                  nFactors;                           // Total num of factors
    double                                  lambdaPrior_reg;                    // Weight for regularization prior -> neighbour cells have similar wind vectors
    double                                  lambdaPrior_mass_conservation;      // Weight for mass conservation law prior
    double                                  lambdaPrior_obstacles;              // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind

    struct TobservationGMRF
    {
        size_t  cell_idx;
        double	windX;
        double	windY;
        double	lambda;
        bool	time_invariant;						//if the observation will lose weight (lambda) as time goes on (default false)
    };

    //GMRF structures
    std::vector<Eigen::Triplet<double> >        J;              // the Jacobian
    std::vector<Eigen::Triplet<double> >        Lambda;         // the information matrix (weights)
    std::vector<TobservationGMRF>               activeObs;      // Vector with the active observations and their respective Information

    //functions
    bool        is_cell_free(size_t id_gmrf);
    bool        check_connectivity_between2cells(size_t idx_1_gmrf, size_t idx_2_gmrf);
    inline int  x2idx(float x) const { return static_cast<int>( (x-m_x_min)/m_resolution ); }
    inline int  y2idx(float y) const { return static_cast<int>( (y-m_y_min)/m_resolution ); }
    inline int  xy2idx(float x,float y) const { return x2idx(x) + y2idx(y)*m_size_x; }
    void        id2cellxy(size_t id, size_t &cell_x, size_t &cell_y);
    void        id2xy(size_t id, double &x, double &y);

    //Visualization    
    void save_grmf_factor_graph(std::vector<Eigen::Triplet<double> > &Jout, std::vector<Eigen::Triplet<double> > &Aout, Eigen::VectorXd &yout);
    void save_grmf_factor_graph(Eigen::SparseMatrix<double> &H, Eigen::VectorXd &G);
    visualization_msgs::msg::Marker line_list, line_list_obs;
    void init_colormaps(std::string colormap);
    float color_r[200];
    float color_g[200];
    float color_b[200];

    //pcl::PointCloud<pcl::PointXYZRGB> template_cells[NUM_CELL_TEMPLATES];
    //void init_pcl_templates(std::string colormap, int max_points_cell);

    //float                                   GMRF_lambdaPrior;		//!< The information (Lambda) of fixed map constraints
    //std::vector<Eigen::Triplet<double> >    H_prior;        // the prior part of H
    //Eigen::VectorXd                         g;              // Gradient vector
    //std::multimap<size_t,size_t>            cell_interconnections;		//Store the interconnections (relations) of each cell with its neighbourds
};
