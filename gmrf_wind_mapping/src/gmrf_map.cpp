#include "gmrf_map.h"




/*---------------------------------------------------------------
                        Constructor
  ---------------------------------------------------------------*/
CGMRF_map::CGMRF_map(const nav_msgs::OccupancyGrid &oc_map,
                     float cell_size,
                     double m_lambdaPrior_reg,
                     double m_lambdaPrior_mass_conservation,
                     double m_lambdaPrior_obstacles,
                     std::string m_colormap, int max_points_cell)
{
    try
    {
        // Copy params to internal variables
        m_Ocgridmap = oc_map;                       // Occupancy gridMap from MapServer
        m_resolution = cell_size;                   // Desired resolution to build the GMRF (m)
        lambdaPrior_reg = m_lambdaPrior_reg;
        lambdaPrior_mass_conservation = m_lambdaPrior_mass_conservation;
        lambdaPrior_obstacles = m_lambdaPrior_obstacles;

        //Set initial GMRF dimensions as the OccupancyMap (in meters)
        double x_min =  oc_map.info.origin.position.x;
        double x_max =  oc_map.info.origin.position.x + oc_map.info.width*oc_map.info.resolution;
        double y_min =  oc_map.info.origin.position.y;
        double y_max =  oc_map.info.origin.position.y + oc_map.info.height*oc_map.info.resolution;

        // Adjust size to complaint with the desired resolution (m_resolution):
        m_x_min = m_resolution*round(x_min/m_resolution);
        m_y_min = m_resolution*round(y_min/m_resolution);
        m_x_max = m_resolution*round(x_max/m_resolution);
        m_y_max = m_resolution*round(y_max/m_resolution);


        // Now the number of cells should be integers:
        m_size_x = round((m_x_max-m_x_min)/m_resolution);
        m_size_y = round((m_y_max-m_y_min)/m_resolution);
        N = m_size_x*m_size_y;


        //For visualization only - pot lines between connected nodes
        line_list.header.stamp = ros::Time();
        line_list.ns = "factors_reg";
        line_list.type = visualization_msgs::Marker::LINE_LIST;
        line_list.action = visualization_msgs::Marker::ADD;
        line_list.id = 0;
        line_list.points.clear();
        //color = blue
        line_list.scale.x = 0.02;
        line_list.color.r = 0.0;
        line_list.color.g = 0.0;
        line_list.color.b = 1.0;
        line_list.color.a = 1.0;

        //For visualization only - plot lines at nodes with Wx==0 or Wy==0 (factors due to obstacles)
        line_list_obs.header.stamp = ros::Time();
        line_list_obs.ns = "factors_obs";
        line_list_obs.type = visualization_msgs::Marker::LINE_LIST;
        line_list_obs.action = visualization_msgs::Marker::ADD;
        line_list_obs.id = 1;
        line_list_obs.points.clear();
        //color = blue
        line_list_obs.scale.x = 0.02;
        line_list_obs.color.r = 1.0;
        line_list_obs.color.g = 0.0;
        line_list_obs.color.b = 0.0;
        line_list_obs.color.a = 1.0;

        //points to add
        geometry_msgs::Point p;
        p.x = 0;
        p.y = 0;
        p.z = 0;


        //-------------------------------------------------------
        // INIT RANDOM FIELD AND CREATE CONNEXIONS BETWEEN NODES
        //-------------------------------------------------------
        ROS_INFO("[CGMRF] Generating GMRF for WIND estimation' ");

        //1. Init the map container
        //-------------------------
        TRandomFieldCell init_cell;
        init_cell.mean = 0.0;
        init_cell.std = 0.0;
        m_map.assign(2*N, init_cell);   //Since we have Wx and Wy, we refer to them as: Wx in the range [0,N-1], Wy in the range [N,2N-1]

        ROS_INFO("--------------------------------");
        ROS_INFO("[CGMRF] GMRF created:");
        ROS_INFO("Occupancy size: x=(%.2f,%.2f)[m] y=(%.2f,%.2f)[m]", x_min, x_max, y_min, y_max);
        ROS_INFO("Occupancy size: (%u,%u) cells with cell_size %.2fm", oc_map.info.width, oc_map.info.height, oc_map.info.resolution);
        ROS_INFO("GMRF size:      x=(%.2f,%.2f)[m] y=(%.2f,%.2f)[m]", m_x_min, m_x_max, m_y_min, m_y_max);
        ROS_INFO("GMRF size:      (%lu,%lu) cells with cell_size %.2fm", m_size_x, m_size_y, m_resolution);
        ROS_INFO("GMRF size:      N = %lu cells, 2N = %lu nodes", N, m_map.size());
        ROS_INFO("--------------------------------");


        //2. Set NumFactors
        //-------------------------
        // Determine number of "static" connections=factors between the nodes in the GMRF
        // We use the occupancy_map to that end, with the intention to allocate memory in the sparse matrices(speed)
        nPriorFactors = 2*( (m_size_x -1) * m_size_y + m_size_x * (m_size_y -1) );	//Approximation: L = (Nr-1)*Nc + Nr*(Nc-1); Full connected

        /* Proper way....but slower?
        //2.1 Num prior factors (L)
        nPriorFactors = 0;
        for (size_t j=0; j<N; j++)
        {
            //Get (r,c) %(row, col)
            jc = mod(j,Nc);
            if (jc == 0)
                jc=Nc;
            jr = ceil(j/Nc);


            // Factor with the right node
            //-----------------------------
            if (jc<Nc)
            {
                // if Cj and Cj+1 are free
                if( metricmap(j)==0 && metricmap(j+1)==0 )
                {
                    if (factor_regularization)
                        nPriorFactors = nPriorFactors +2;   //for Wx and Wy
                }
                // if Cj or Cj+1 are occupied --> Factor of Null perpendicular value
                else if( metricmap(j)==0 || metricmap(j+1)==0 )
                {
                    if (factor_tangencial_to_obstacles)
                        nPriorFactors = nPriorFactors +1;
                // if both are occupied --> No factor
                }
            }


            // Factor with the above node
            //---------------------------
            if (jr<Nr)
                // if Cj and Cj+Nc are free
                if( metricmap(j)==0 && metricmap(j+Nc)==0 )
                     if factor_regularization
                        nPriorFactors = nPriorFactors +2;   //for Wx and Wy
                    end;

                // if Cj or Cj+Nc are occupied --> Factor of Null perpendicular value
                elseif( metricmap(j)==0 || metricmap(j+Nc)==0 )
                    if factor_tangencial_to_obstacles
                        nPriorFactors = nPriorFactors +1;
                    end;

                // if both are occupied --> No factor
                end;
            end;


            if (factor_mass_conservation_law)
                // Only set the factor if cell "j" is free, and it is not in the contour of the map
                if (metricmap(j)==0 && jc>1 && jc<Nc && jr>1 && jr<Nr)
                    %As soon as any of its neighbour cells is free, set the factor
                    if (metricmap(j-1)==0)
                        nPriorFactors = nPriorFactors + 1;
                    elseif (metricmap(j+1)==0)
                        nPriorFactors = nPriorFactors + 1;
                    elseif (metricmap(j-Nc)==0 )
                        nPriorFactors = nPriorFactors + 1;
                    elseif (metricmap(j+Nc)==0)
                        nPriorFactors = nPriorFactors + 1;
                    end;
                end;
            end;
        }
        */

        //2.2 Num Observation factors (M)
        nObsFactors = 0;

        //2.3 TOTAL Number of Factors
        nFactors = nPriorFactors + nObsFactors;
        ROS_INFO("[CGMRF] Setting up Prior-factors");



        //3. Reserve memory to SpeedUp
        //-----------------------------
        J.clear();
        J.reserve(5*nFactors);          //Each factor accounts for the connectivity of multiple nodes -->multiple entries
        Lambda.clear();
        Lambda.reserve(nFactors);       //Diagonal -> 1 entry for each factor



        //------------------------------------
        //4. Build Prior Factors (just once)
        //-------------------------------------
        //Given the possibility of using different cell_sizes for the GMRF and Occupancy maps,
        //We need to employ the Occupancy map to determine the real interconnections between cells and therefore create the factors between nodes in the GMRF.

        size_t count = 0;
        for (size_t j=0; j<N; j++)		//For each cell in the GMRF
        {
            // Get cell_x and cell_y in the GMRF representation
            size_t jx,jy;
            id2cellxy(j, jx, jy);

            if (!is_cell_free(j))
            {
                // Force occupied cell j to 0 value
                // This is mandatory for the correct solution of the system
                Eigen::Triplet<double> lambda_entry(count,count, lambdaPrior_obstacles);
                Eigen::Triplet<double> J_entry(count,j, 1);
                Lambda.push_back(lambda_entry);
                J.push_back(J_entry);
                count++;
                Eigen::Triplet<double> lambda_entry2(count,count, lambdaPrior_obstacles);
                Eigen::Triplet<double> J_entry2(count,j+N, 1);
                Lambda.push_back(lambda_entry2);
                J.push_back(J_entry2);
                count++;
            }

            //Factor with the right node: (j <--> j+1)
            //-----------------------------------------
            if (jx < (m_size_x-1))
            {
                if( is_cell_free(j) && is_cell_free(size_t(j+1)) )
                {
                    if (check_connectivity_between2cells(j,size_t(j+1)))
                    {
                        //Create a regularization factor to link both nodes
                        // Wx range [1,N]
                        Eigen::Triplet<double> lambda_entry(count,count, lambdaPrior_reg);
                        Eigen::Triplet<double> J_entry1(count,j, 1);
                        Eigen::Triplet<double> J_entry2(count,size_t(j+1), -1);
                        Lambda.push_back(lambda_entry);
                        J.push_back(J_entry1);
                        J.push_back(J_entry2);
                        count++;

                        // Wy range [N+1,2N]
                        Eigen::Triplet<double> lambda_entry2(count,count, lambdaPrior_reg);
                        Eigen::Triplet<double> J_entry3(count,j+N, 1);
                        Eigen::Triplet<double> J_entry4(count,size_t(j+N+1), -1);
                        Lambda.push_back(lambda_entry2);
                        J.push_back(J_entry3);
                        J.push_back(J_entry4);
                        count++;

                        //plot marker between cell j and j+1
                        id2xy(j, p.x, p.y);
                        p.z = 0;
                        line_list.points.push_back(p);
                        id2xy(j+1, p.x, p.y);
                        line_list.points.push_back(p);
                    }
                    else
                    {
                        // An obstacle is in between both cells
                        //1. Do not create a regularization link
                        //2. Force Wx=0 at both cells
                        Eigen::Triplet<double> lambda_entry(count,count, lambdaPrior_obstacles);
                        Eigen::Triplet<double> J_entry(count,j, 1);
                        Lambda.push_back(lambda_entry);
                        J.push_back(J_entry);
                        count++;
                        //Plot marker to display Wx=0
                        id2xy(j, p.x, p.y);
                        p.z = 0;
                        p.x += m_resolution/5;
                        p.y += m_resolution/5;
                        line_list_obs.points.push_back(p);
                        p.y -= 2* m_resolution/5;
                        line_list_obs.points.push_back(p);

                        Eigen::Triplet<double> lambda_entry2(count,count, lambdaPrior_obstacles);
                        Eigen::Triplet<double> J_entry2(count,j+1, 1);
                        Lambda.push_back(lambda_entry2);
                        J.push_back(J_entry2);
                        count++;
                        //Plot marker to display Wx=0
                        id2xy(j+1, p.x, p.y);
                        p.z = 0;
                        p.x -= m_resolution/5;
                        p.y += m_resolution/5;
                        line_list_obs.points.push_back(p);
                        p.y -= 2* m_resolution/5;
                        line_list_obs.points.push_back(p);
                    }
                }
                else if (is_cell_free(j))
                {
                    // Force Wx=0 at j
                    Eigen::Triplet<double> lambda_entry(count,count, lambdaPrior_obstacles);
                    Eigen::Triplet<double> J_entry(count,j, 1);
                    Lambda.push_back(lambda_entry);
                    J.push_back(J_entry);
                    count++;
                    //Plot marker to display Wx=0
                    id2xy(j, p.x, p.y);
                    p.z = 0;
                    p.x += m_resolution/5;
                    p.y += m_resolution/5;
                    line_list_obs.points.push_back(p);
                    p.y -= 2* m_resolution/5;
                    line_list_obs.points.push_back(p);
                }
                else if (is_cell_free(j+1))
                {
                    // Force Wx=0 at j+1
                    Eigen::Triplet<double> lambda_entry(count,count, lambdaPrior_obstacles);
                    Eigen::Triplet<double> J_entry(count,j+1, 1);
                    Lambda.push_back(lambda_entry);
                    J.push_back(J_entry);
                    count++;
                    //Plot marker to display Wx=0
                    id2xy(j+1, p.x, p.y);
                    p.z = 0;
                    p.x -= m_resolution/5;
                    p.y += m_resolution/5;
                    line_list_obs.points.push_back(p);
                    p.y -= 2* m_resolution/5;
                    line_list_obs.points.push_back(p);
                }
                //else --> Both cells occupied -> Do nothing!
            }


            //Factor with the upper node: (j <--> j+m_size_x)
            //------------------------------------------------
            if (jy < (m_size_y-1))
            {
                if( is_cell_free(j) && is_cell_free(j+m_size_x) )
                {
                    if (check_connectivity_between2cells(j,j+m_size_x))
                    {
                        //Create a regularization factor to link both nodes
                        // Wx range [1,N]
                        Eigen::Triplet<double> lambda_entry(count,count, lambdaPrior_reg);
                        Eigen::Triplet<double> J_entry1(count,j, 1);
                        Eigen::Triplet<double> J_entry2(count,j+m_size_x, -1);
                        Lambda.push_back(lambda_entry);
                        J.push_back(J_entry1);
                        J.push_back(J_entry2);
                        count++;

                        // Wy range [N+1,2N]
                        Eigen::Triplet<double> lambda_entry2(count,count, lambdaPrior_reg);
                        Eigen::Triplet<double> J_entry3(count,j+N, 1);
                        Eigen::Triplet<double> J_entry4(count,j+N+m_size_x, -1);
                        Lambda.push_back(lambda_entry2);
                        J.push_back(J_entry3);
                        J.push_back(J_entry4);
                        count++;

                        //plot marker between cell j and j+m_size_x
                        id2xy(j, p.x, p.y);
                        p.z = 0;
                        line_list.points.push_back(p);
                        id2xy(j+m_size_x, p.x, p.y);
                        line_list.points.push_back(p);
                    }
                    else
                    {
                        // An obstacle is in between both cells
                        //1. Do not create a regularization link
                        //2. Force Wy=0 at both cells
                        // Force Wy=0
                        Eigen::Triplet<double> lambda_entry(count,count, lambdaPrior_obstacles);
                        Eigen::Triplet<double> J_entry(count,j+N, 1);
                        Lambda.push_back(lambda_entry);
                        J.push_back(J_entry);
                        count++;
                        //Plot marker to display Wy=0
                        id2xy(j, p.x, p.y);
                        p.z = 0;
                        p.y += m_resolution/5;
                        p.x -= m_resolution/5;
                        line_list_obs.points.push_back(p);
                        p.x += 2* m_resolution/5;
                        line_list_obs.points.push_back(p);

                        // Force Wy=0
                        Eigen::Triplet<double> lambda_entry2(count,count, lambdaPrior_obstacles);
                        Eigen::Triplet<double> J_entry2(count,j+N+m_size_x, 1);
                        Lambda.push_back(lambda_entry2);
                        J.push_back(J_entry2);
                        count++;
                        //Plot marker to display Wy=0
                        id2xy(j+m_size_x, p.x, p.y);
                        p.z = 0;
                        p.y -= m_resolution/5;
                        p.x -= m_resolution/5;
                        line_list_obs.points.push_back(p);
                        p.x += 2* m_resolution/5;
                        line_list_obs.points.push_back(p);
                    }
                }
                else if (is_cell_free(j))
                {
                    // Force Wy=0 at j
                    Eigen::Triplet<double> lambda_entry(count,count, lambdaPrior_obstacles);
                    Eigen::Triplet<double> J_entry(count,j+N, 1);
                    Lambda.push_back(lambda_entry);
                    J.push_back(J_entry);
                    count++;
                    //Plot marker to display Wy=0
                    id2xy(j, p.x, p.y);
                    p.z = 0;
                    p.y += m_resolution/5;
                    p.x -= m_resolution/5;
                    line_list_obs.points.push_back(p);
                    p.x += 2* m_resolution/5;
                    line_list_obs.points.push_back(p);
                }
                else if (is_cell_free(j+m_size_x))
                {
                    // Force Wy=0 at j+m_size_x
                    Eigen::Triplet<double> lambda_entry2(count,count, lambdaPrior_obstacles);
                    Eigen::Triplet<double> J_entry2(count,j+N+m_size_x, 1);
                    Lambda.push_back(lambda_entry2);
                    J.push_back(J_entry2);
                    count++;
                    //Plot marker to display Wy=0
                    id2xy(j+m_size_x, p.x, p.y);
                    p.z = 0;
                    p.y -= m_resolution/5;
                    p.x -= m_resolution/5;
                    line_list_obs.points.push_back(p);
                    p.x += 2* m_resolution/5;
                    line_list_obs.points.push_back(p);
                }
                //else --> Both cells occupied -> Do nothing!
            }


            // Factors for mass conservative Law (avoid borders of map)
            //---------------------------------------------------------

            if (is_cell_free(j) && jx>0 && jx<m_size_x-1 && jy>0 && jy<m_size_y-1)
            {
                //As soon as any of its 8 clossest neighbour cells is free, set the factor
                bool set = false;
                if (is_cell_free(j-1)){
                    Eigen::Triplet<double> J_entry(count,j-1, -1);
                    J.push_back(J_entry);
                    set = true;
                }
                if (is_cell_free(j+1)){
                    Eigen::Triplet<double> J_entry(count,j+1, 1);
                    J.push_back(J_entry);
                    set = true;
                }
                if (is_cell_free(j-m_size_x)){
                    Eigen::Triplet<double> J_entry(count,j+N-m_size_x, -1);
                    J.push_back(J_entry);
                    set = true;
                }
                if (is_cell_free(j+m_size_x)){
                    Eigen::Triplet<double> J_entry(count,j+N+m_size_x, 1);
                    J.push_back(J_entry);
                    set = true;
                }

                //Diagonals
                if (is_cell_free(j+m_size_x-1)){
                    Eigen::Triplet<double> J_entry(count,j+m_size_x-1, -0.5);
                    Eigen::Triplet<double> J_entry2(count,j+m_size_x-1+N, 0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }
                if (is_cell_free(j+m_size_x+1)){
                    Eigen::Triplet<double> J_entry(count,j+m_size_x+1, 0.5);
                    Eigen::Triplet<double> J_entry2(count,j+m_size_x+1+N, 0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }
                if (is_cell_free(j-m_size_x+1)){
                    Eigen::Triplet<double> J_entry(count,j-m_size_x+1, 0.5);
                    Eigen::Triplet<double> J_entry2(count,j-m_size_x+1+N, -0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }
                if (is_cell_free(j-m_size_x-1)){
                    Eigen::Triplet<double> J_entry(count,j-m_size_x-1, -0.5);
                    Eigen::Triplet<double> J_entry2(count,j-m_size_x-1+N, -0.5);
                    J.push_back(J_entry);
                    J.push_back(J_entry2);
                    set = true;
                }

                //If factor is to be set...
                if (set){
                    Eigen::Triplet<double> lambda_entry(count,count, lambdaPrior_mass_conservation);
                    Lambda.push_back(lambda_entry);
                    count++;
                    set = false;
                }
            }

        } // end for setting factors


        //Set number of Factors
        nPriorFactors = count;
        nFactors = nPriorFactors + nObsFactors;
        activeObs.clear();
        ROS_INFO("[CGMRF] Initialization complete: %lu factors for a map size of 2N=%lu nodes", nFactors, m_map.size());

        //Set the colormap to display the wind vectors
        init_colormaps("jet");


        // DEBUG: Save to file
        //-------------------
        /*
        Eigen::VectorXd y_empty;
        y_empty.resize(nFactors);
        y_empty.fill(0.0);
        save_grmf_factor_graph(J,Lambda,y_empty);
        */


        //DEGUB : ADD FIXED OBSERVATION
        /*
        TobservationGMRF new_obs;
        const int cellIdx = xy2idx( 3.0, 3.0 );
        new_obs.cell_idx = cellIdx;
        new_obs.windX = 0.10;
        new_obs.windY = 1.0;
        new_obs.lambda = 13;
        new_obs.time_invariant = false;		//Default behaviour, the obs will lose weight with time.
        ROS_INFO("[GMRF] DEMO obs: Wx = %.2f m/s Wy = %.2f m/s at cell %lu\n\n", new_obs.windX,new_obs.windY,new_obs.cell_idx);
        activeObs.push_back(new_obs);
        nObsFactors += 2;    //we add 2 factors for each observation to account for Wx and Wy components
        */


    }catch(std::exception e){
        ROS_ERROR("[GMRF-Constructor] Exception: %s ", e.what() );
    }
}


/*---------------------------------------------------------------
                        Destructor
  ---------------------------------------------------------------*/
CGMRF_map::~CGMRF_map(){}



/*---------------------------------------------------------------
                        Cell index transformations
  ---------------------------------------------------------------*/
// Get x,y in cells (in the GMRF representation) from the general index in the array
void CGMRF_map::id2cellxy(size_t id, size_t &cell_x, size_t &cell_y)
{
    cell_x = id % m_size_x;
    cell_y = (size_t) floor(id/m_size_x);
}

// Get pose x,y (in meters) (in the GMRF representation) from the general index in the array
void CGMRF_map::id2xy(size_t id, double &x, double &y)
{
    size_t cell_x, cell_y;
    id2cellxy(id,cell_x,cell_y);

    x = m_x_min + (cell_x*m_resolution) + (m_resolution/2);
    y = m_y_min + (cell_y*m_resolution) + (m_resolution/2);
}


/*---------------------------------------------------------------
             Check if a cell is free of obstacles
  ---------------------------------------------------------------*/
//Check at OccupancyMap level, if a cell is free of obstacles (checking the cell center at GMRF resolution)
bool CGMRF_map::is_cell_free(size_t id_gmrf)
{
    //The pose x,y (meters) of cell center
    double cell_1_x,cell_1_y;
    id2xy(id_gmrf,cell_1_x,cell_1_y);

    //Get corresponding cell_idx in the Occupancy Gridmap
    //IMPORTANT --> Use the resolution and size of Occupancy Gridmap (not the GMRF)
    int id_oc;
    id_oc = static_cast<int>( (cell_1_x-m_Ocgridmap.info.origin.position.x)/m_Ocgridmap.info.resolution );  //x component
    id_oc += static_cast<int>( (cell_1_y-m_Ocgridmap.info.origin.position.y)/m_Ocgridmap.info.resolution ) * m_Ocgridmap.info.width;    //y component


    /* DEBUG
    double cellx,celly;
    cellx = m_Ocgridmap.info.origin.position.x + (id_oc % m_Ocgridmap.info.width)*m_Ocgridmap.info.resolution + m_Ocgridmap.info.resolution/2;
    celly = m_Ocgridmap.info.origin.position.y + ((size_t) floor(id_oc/m_Ocgridmap.info.width))*m_Ocgridmap.info.resolution + m_Ocgridmap.info.resolution/2;
    */

    try
    {
        //Check occupancy
        if (m_Ocgridmap.data[id_oc]>=50.0)
        {
            //ROS_INFO("[GMRF] OCCUPIED %lu = (%.2f,%.2f) --> %lu in Occ",idx_1_gmrf,cell_1_x,cell_1_y, idx_1_oc);
            return false;
        }
        else
        {
            //ROS_INFO("[GMRF] FREE %lu = (%.2f,%.2f) --> %lu in Occ",idx_1_gmrf,cell_1_x,cell_1_y, idx_1_oc);
            return true;
        }
    }
    catch(std::exception e){
        ROS_ERROR("[GMRF] Exception while checking cell freedom: %s ", e.what() );
    }
}


/*---------------------------------------------------------------
             Check cell interconnectivity
  ---------------------------------------------------------------*/
//Check at OccupancyMap level, if two cells are interconnected, that is no obstacles in between them.
//If ture, we will set a regularization factor.
bool CGMRF_map::check_connectivity_between2cells(size_t idx_1_gmrf, size_t idx_2_gmrf)
{
    try
    {
        //Ge poses (x,y) of the cell centers in GMRF map
        double cell_1_x,cell_1_y,cell_2_x,cell_2_y;
        id2xy(idx_1_gmrf,cell_1_x,cell_1_y);
        id2xy(idx_2_gmrf,cell_2_x,cell_2_y);

        //Get corresponding cell_idx in the Occupancy Gridmap
        //IMPORTANT --> Use the resolution and size of Occupancy Gridmap (not the GMRF)
        int idx_1_oc, idx_2_oc;
        idx_1_oc = static_cast<int>( (cell_1_x-m_Ocgridmap.info.origin.position.x)/m_Ocgridmap.info.resolution );  //x component
        idx_1_oc += static_cast<int>( (cell_1_y-m_Ocgridmap.info.origin.position.y)/m_Ocgridmap.info.resolution ) * m_Ocgridmap.info.width;    //y component
        idx_2_oc = static_cast<int>( (cell_2_x-m_Ocgridmap.info.origin.position.x)/m_Ocgridmap.info.resolution );  //x component
        idx_2_oc += static_cast<int>( (cell_2_y-m_Ocgridmap.info.origin.position.y)/m_Ocgridmap.info.resolution ) * m_Ocgridmap.info.width;    //y component

        //check if cells are in the same row of the GMRF map
        bool horizontal = false;
        if (idx_2_gmrf == idx_1_gmrf+1)
            horizontal = true;


        //Check that a straigh line between both cells centers is free of obstacles
        bool connected = true;
        if (horizontal)
        {
            for (size_t p=idx_1_oc; p<idx_2_oc; p++)
            {
                if (m_Ocgridmap.data[p]>=60.0)
                {
                    connected = false;
                    break;
                }
            }
        }
        else
        {
            for (size_t p=idx_1_oc; p<idx_2_oc; p+=m_Ocgridmap.info.width)
            {
                if (m_Ocgridmap.data[p]>=50.0)
                {
                    connected = false;
                    break;
                }
            }
        }

        return connected;
    }
    catch(std::exception e){
        ROS_ERROR("[GMRF] Exception while checking cells interconnections: %s ", e.what() );
    }
}




/*---------------------------------------------------------------
             Insert New Wind Observation
---------------------------------------------------------------*/
void CGMRF_map::insertObservation_GMRF(double wind_speed, double wind_direction, double x_pos, double y_pos, double lambdaObs)
{
    try
    {
        // Fill new Observation
        // The wind vector provided is already the DownWind direction in the /map reference system
        const int cellIdx = xy2idx( x_pos, y_pos );
        TobservationGMRF new_obs;
        new_obs.cell_idx = cellIdx;
        new_obs.windX = wind_speed * cos(wind_direction);
        new_obs.windY = wind_speed * sin(wind_direction);
        new_obs.lambda = lambdaObs;
        new_obs.time_invariant = false;		//Default behaviour, the obs will lose weight with time.
        ROS_INFO("[GMRF] New obs: Wx = %.2f m/s Wy = %.2f m/s", new_obs.windX,new_obs.windY);

        // Add Observation to GMRF
        activeObs.push_back(new_obs);
        nObsFactors += 2;    //we add 2 factors foe each observation to account for Wx and Wy components

        // NOTE --> We create 4 observations to expand a bit the measurement impact, replicating the content to neighbour cells
        if (is_cell_free(cellIdx-1))
        {
            new_obs.cell_idx = cellIdx-1;
            activeObs.push_back(new_obs);
            nObsFactors += 2;
        }
        if (is_cell_free(cellIdx-m_size_x))
        {
            new_obs.cell_idx = cellIdx-m_size_x;
            activeObs.push_back(new_obs);
            nObsFactors += 2;
        }
        if (is_cell_free(cellIdx-m_size_x-1))
        {
            new_obs.cell_idx = cellIdx-m_size_x-1;
            activeObs.push_back(new_obs);
            nObsFactors += 2;
        }

    }catch(std::exception e){
        ROS_ERROR("[GMRF] Exception while Inserting new Observation: %s ", e.what() );
    }
}




/*---------------------------------------------------------------
                    updateMapEstimation_GMRF
  ---------------------------------------------------------------*/
void  CGMRF_map::updateMapEstimation_GMRF(float lambdaObsLoss)
{
    try
    {
        /*
         * J (Jacobian) The J matrix contains the dr/dm for every factor in the graph
         *              J is size (nFactors x NumNodes)
         *
         * Lambda (weights) Is the Diagonal information matrix (contains the weights for each factor)
         *              Lambda is size (nFactors x nFactors)
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

        //1. Get current number of factors (nPriorFactors is constant, but nObsFactors is dynamic)
        nFactors = nPriorFactors + nObsFactors;

        //2. Copy The prior part of Jacobian and Lambda matrices
        std::vector<Eigen::Triplet<double> > J_temp;
        J_temp.reserve( J.size()+nObsFactors );
        std::copy( J.begin(), J.end(), back_inserter(J_temp) );

        std::vector<Eigen::Triplet<double> > Lambda_temp;
        Lambda_temp.reserve( Lambda.size()+nObsFactors );
        std::copy( Lambda.begin(), Lambda.end(), back_inserter(Lambda_temp) );

        Eigen::VectorXd y_temp;
        y_temp.resize(nFactors);
        y_temp.fill(0.0);



        //3. Include Active Observations into Jacobian and Lambda
        size_t count = nPriorFactors;                                   // start after the already introduced prior factors
        for (std::vector<TobservationGMRF>::iterator ito = activeObs.begin(); ito !=activeObs.end(); ++ito)
        {
            //Each observation translates to 2 factors (Wx,Wy)
            // Wx range [1,N]
            Eigen::Triplet<double> lambda_entry(count,count, ito->lambda);
            Eigen::Triplet<double> J_entry(count,ito->cell_idx, 1);
            Lambda_temp.push_back(lambda_entry);
            J_temp.push_back(J_entry);
            y_temp[count] = ito->windX;
            count++;

            // Wy range [N+1,2N]
            Eigen::Triplet<double> lambda_entry2(count,count, ito->lambda);
            Eigen::Triplet<double> J_entry2(count,ito->cell_idx+N, 1);
            Lambda_temp.push_back(lambda_entry2);
            J_temp.push_back(J_entry2);
            y_temp[count] = ito->windY;
            count++;
        }

        // DEBUG - Save to file
        //save_grmf_factor_graph(J_temp,Lambda_temp,y_temp);


        //3. Build Matrices (J, J', A, H, G)
        Eigen::SparseMatrix<double> Jsparse(nFactors,2*N);				// declares a column-major sparse matrix type of float
        Jsparse.setFromTriplets(J_temp.begin(), J_temp.end() );
        ROS_INFO("          [GMRF] Jsparse is (%u,%u)",Jsparse.rows(),Jsparse.cols());

        Eigen::SparseMatrix<double> JsparseT;//(2*N,nFactors);				// declares a column-major sparse matrix type of float
        JsparseT = Eigen::SparseMatrix<double>(Jsparse.transpose());
        ROS_INFO("          [GMRF] JsparseT is (%u,%u)",JsparseT.rows(),JsparseT.cols());

        Eigen::SparseMatrix<double> Asparse(nFactors,nFactors);			// declares a column-major sparse matrix type of float
        Asparse.setFromTriplets(Lambda_temp.begin(), Lambda_temp.end() );
        ROS_INFO("          [GMRF] Asparse is (%u,%u)",Asparse.rows(),Asparse.cols());

        Eigen::SparseMatrix<double> Hsparse;//(2*N,2*N);   				// declares a column-major sparse matrix type of float
        Hsparse = JsparseT * Asparse * Jsparse;
        ROS_INFO("          [GMRF] Hsparse is (%u,%u)",Hsparse.rows(),Hsparse.cols());

        Eigen::VectorXd G = JsparseT * Asparse * y_temp;
        ROS_INFO("          [GMRF] G is (%lu,%lu)",G.rows(),G.cols());


        // DEBUG - Save to file
        //save_grmf_factor_graph(Hsparse, G);

        //----------
        // 4. SOLVE
        //----------
        // We need to solve: H * inc_m = -G
        // In an interative scenario: m = m + inc_m;
        // In our case, we do not need to consider the previous state, so m = inc_m
        // We use a Cholesky Factorization of Hessian --> chol( P * H * inv(P) )
        Eigen::SimplicialLLT< Eigen::SparseMatrix <double> > solver;
        solver.compute(Hsparse);
        Eigen::VectorXd m_inc = solver.solve(G);

        ROS_INFO("[GMRF] system solved with solution size (%lu,%lu)", m_inc.rows(), m_inc.cols());
        std::ofstream file("/home/jgmonroy/gmrf_solution.txt");
        if (file.is_open())
        {
            file << m_inc;
        }
        file.close();

        // 5. Update GMRF values from current solution
        for (size_t j=0; j<m_map.size(); j++)
        {
            m_map[j].mean = m_inc(j);       // Not iterative! no need to increment previous state
            m_map[j].std = 0.0;             // Not used right now
        }



        // 7. Update Information/Strength of Active Observations
        //------------------------------------------------------- 
        std::vector<TobservationGMRF>::iterator ito = activeObs.begin();
        while ( ito!=activeObs.end() )
        {
            if (ito->time_invariant==false)
            {
                ito->lambda -= lambdaObsLoss;
                if (ito->lambda <= 0.0)
                {
                    ito = activeObs.erase(ito);
                    nObsFactors -= 2;
                }
                else
                    ++ito;
            }else
                ++ito;
        }
        ROS_INFO("[GMRF] %lu ObservationFactors are active", nObsFactors);
    }catch(std::exception e){
        ROS_ERROR("[GMRF] Exception Updating the maps: %s ", e.what() );
    }
}


Eigen::Vector2d CGMRF_map::getEstimation(double x, double y){
    int i = xy2idx(x,y);
    double module = sqrt(pow(m_map[i].mean,2) + pow(m_map[i+N].mean,2));
    double direction = atan2(m_map[i+N].mean, m_map[i].mean);

    return Eigen::Vector2d(module, direction);
}

void CGMRF_map::get_as_markerArray(visualization_msgs::MarkerArray &ma, std::string frame_id)
{
    ma.markers.clear();
    //options (debug)
    bool plot_cell_centers = false;
    bool plot_factors = false;
    bool plot_wind_vectors = true;


    if( plot_cell_centers)
    {
        //marker-points at all cells
        visualization_msgs::Marker marker_free;
        visualization_msgs::Marker marker_occ;
        marker_free.header.frame_id = marker_occ.header.frame_id = frame_id.c_str();
        marker_free.header.stamp = marker_occ.header.stamp = ros::Time();
        marker_free.ns = marker_occ.ns = "cell_centers";
        marker_free.type = marker_occ.type = visualization_msgs::Marker::POINTS;
        marker_free.action = marker_occ.action = visualization_msgs::Marker::ADD;
        marker_free.id = 2;
        marker_occ.id = 3;
        // POINTS markers use x and y scale for width/height respectively
        marker_free.scale.x = 0.1;
        marker_free.scale.y = 0.1;
        marker_free.color.r = 0.0;
        marker_free.color.g = 1.0;
        marker_free.color.b = 0.0;
        marker_free.color.a = 1.0;

        marker_occ.scale.x = 0.1;
        marker_occ.scale.y = 0.1;
        marker_occ.color.r = 1.0;
        marker_occ.color.g = 0.0;
        marker_occ.color.b = 0.0;
        marker_occ.color.a = 1.0;

        //fill points
        marker_free.points.clear();
        marker_occ.points.clear();
        for (size_t i=0; i<N; i++)
        {
            double cell_center_x, cell_center_y;
            geometry_msgs::Point p;
            id2xy(i, cell_center_x, cell_center_y);
            p.x = cell_center_x;
            p.y = cell_center_y;

            if (is_cell_free(i))
                marker_free.points.push_back(p);
            else
                marker_occ.points.push_back(p);
        }
        //Push PointMarker to array
        ma.markers.push_back(marker_occ);
        ma.markers.push_back(marker_free);
    }


    if( plot_factors)
    {
        //Push Line_list marker to array
        line_list.header.frame_id = frame_id.c_str();
        line_list_obs.header.frame_id = frame_id.c_str();
        ma.markers.push_back(line_list);
        ma.markers.push_back(line_list_obs);
    }


    if (plot_wind_vectors)
    {
        //Add an ARROW marker for each node
        visualization_msgs::Marker marker;
        marker.header.frame_id = frame_id.c_str();
        marker.header.stamp = ros::Time();
        marker.ns = "WindVector";
        marker.type = visualization_msgs::Marker::ARROW;
        marker.action = visualization_msgs::Marker::ADD;

        //Get max wind vector in the map (to normalize the plot)
        double max_module = 0.0;
        for (size_t i=0; i<N; i++)
        {
            if (sqrt(pow(m_map[i].mean,2) + pow(m_map[i+N].mean,2) > max_module) )
                max_module = sqrt(pow(m_map[i].mean,2) + pow(m_map[i+N].mean,2));
        }

        for (size_t i=0; i<N; i++)
        {
            //if (is_cell_free(i))
            {
                double module = sqrt(pow(m_map[i].mean,2) + pow(m_map[i+N].mean,2));
                //ROS_INFO("[GMRF] wind(%lu)=(%.2f,%.2f)m/s",i,m_map[i].mean,m_map[i+N].mean );
                if ( module> 0.001)
                {
                    // Set the pose of the marker.
                    marker.id = i+10;
                    double cell_center_x, cell_center_y;
                    id2xy(i, cell_center_x, cell_center_y);
                    marker.pose.position.x = cell_center_x;
                    marker.pose.position.y = cell_center_y;
                    marker.pose.orientation = tf::createQuaternionMsgFromYaw( atan2(m_map[i+N].mean, m_map[i].mean) );
                    //shape
                    marker.scale.x = m_resolution *  (module/max_module);      // arrow length,
                    marker.scale.y = 0.03;           // arrow width
                    marker.scale.z = 0.05;           // arrow height
                    // color -> must normalize to [0-199]
                    size_t idx_color = 199 * (module/max_module);
                    marker.color.r = color_r[idx_color];
                    marker.color.g = color_g[idx_color];
                    marker.color.b = color_b[idx_color];
                    marker.color.a = 1.0;

                    //Push Arrow to array
                    ma.markers.push_back(marker);
                }
            }
        }//end for
    }
}


void CGMRF_map::save_grmf_factor_graph(std::vector<Eigen::Triplet<double> > &Jout, std::vector<Eigen::Triplet<double> > &Aout, Eigen::VectorXd &yout)
{
    bool save_dense = true;
    bool save_sparse = true;
    if (save_dense)
    {
        //1. Jacobian
        Eigen::SparseMatrix<double> Jsparse(nFactors,2*N);				// declares a column-major sparse matrix type of float
        Jsparse.setFromTriplets(Jout.begin(), Jout.end() );
        Eigen::MatrixXd Jdense = Jsparse.toDense();

        // define the format you want, you only need one instance of this...
        const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
        //ROS_INFO("[GMRF] Saving Factor-Graph to file...");
        std::ofstream file("/home/jgmonroy/gmrf_jacobian_dense.txt");
        if (file.is_open())
        {
            file << Jdense.format(CSVFormat) << '\n';
        }
        file.close();


        //2. Information Matrix
        Eigen::SparseMatrix<double> Asparse(nFactors,nFactors);				// declares a column-major sparse matrix type of float
        Asparse.setFromTriplets(Aout.begin(), Aout.end() );
        Eigen::MatrixXd Adense = Asparse.toDense();
        std::ofstream file2("/home/jgmonroy/gmrf_lambda_dense.txt");
        if (file2.is_open())
        {
            file2 << Adense.format(CSVFormat) << '\n';
        }
        file2.close();
    }

    if(save_sparse)
    {
        ROS_INFO("[GMRF] Saving Factor-Graph (list of triplets) to file...");
        ROS_INFO("[GMRF] Jtriplets(%lu,3), Atriplets(%lu,3), numFactors(%lu)",Jout.size(),Aout.size(),yout.rows());

        //1. Jacobian
        std::ofstream file("/home/jgmonroy/gmrf_Jacobian.txt");
        if (file.is_open())
        {
            file << "# Jacobian of the GMRF: row col value" << "\n";
            file << "# Num_cells_x = " << m_size_x << "\n";
            file << "# Num_cells_y = " << m_size_y << "\n";
            file << "# cell_size = " << m_resolution << "\n";
            for (std::vector<Eigen::Triplet<double> >::iterator it = Jout.begin(); it != Jout.end(); it++)
            {
                file << it->row() << " " << it->col()<< " " << it->value() << "\n";
            }
        }
        file.close();

        //2. Information Matrix
        std::ofstream file2("/home/jgmonroy/gmrf_Lambda.txt");
        if (file2.is_open())
        {
            for (std::vector<Eigen::Triplet<double> >::iterator it = Aout.begin(); it != Aout.end(); it++)
            {
                file2 << it->row() << " " << it->col()<< " " << it->value() << "\n";
            }
        }
        file2.close();

        //3. Save vector of observations
        std::ofstream file3("/home/jgmonroy/gmrf_observations.txt");
        if (file3.is_open())
        {
            file3 << yout;
        }
        file3.close();
    }
}


//Save Sparse matrices to file (for debug)
void CGMRF_map::save_grmf_factor_graph(Eigen::SparseMatrix<double> &H, Eigen::VectorXd &G)
{
    // 1. Hessian
    std::ofstream file("/home/jgmonroy/gmrf_hessian.txt");
    if (file.is_open())
    {
        file << "# Hessian of the GMRF: row col value" << "\n";
        for (int k=0; k < H.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(H,k); it; ++it)
            {
                file << it.row() << " "; // row index
                file << it.col() << " "; // col index (here it is equal to k)
                file << it.value() << "\n";
            }
        }
    }
    file.close();

    //Hdense
    Eigen::MatrixXd Hdense = H.toDense();
    // define the format you want, you only need one instance of this...
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file3("/home/jgmonroy/gmrf_hessian_dense.txt");
    if (file3.is_open())
    {
        file3 << Hdense.format(CSVFormat) << '\n';
    }
    file3.close();

    //2. Gradient
    std::ofstream file2("/home/jgmonroy/gmrf_gradient.txt");
    if (file2.is_open())
    {
        file2 << G;
    }
    file2.close();
}


//------------------------------------------
// Build colormaps for visualization
//------------------------------------------
void CGMRF_map::init_colormaps(std::string colormap)
{
    if(colormap.compare("jet")==0)
    {
        float temp_color_r[200]={0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58,0.60,0.62,0.64,0.66,0.68,0.70,0.72,0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56,0.54,0.52,0.50};
        float temp_color_g[200]={0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58,0.60,0.62,0.64,0.66,0.68,0.70,0.72,0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56,0.54,0.52,0.50,0.48,0.46,0.44,0.42,0.40,0.38,0.36,0.34,0.32,0.30,0.28,0.26,0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08,0.06,0.04,0.02,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00};
        float temp_color_b[200]={0.52,0.54,0.56,0.58,0.60,0.62,0.64,0.66,0.68,0.70,0.72,0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56,0.54,0.52,0.50,0.48,0.46,0.44,0.42,0.40,0.38,0.36,0.34,0.32,0.30,0.28,0.26,0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08,0.06,0.04,0.02,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00};
        for (int ix=0;ix<200;ix++) {
            color_r[ix]=temp_color_r[ix];
            color_g[ix]=temp_color_g[ix];
            color_b[ix]=temp_color_b[ix];
        }
    }
    else if(colormap.compare("hot")==0){
        float temp_color_r[200]={0.01,0.03,0.04,0.05,0.07,0.08,0.09,0.11,0.12,0.13,0.15,0.16,0.17,0.19,0.20,0.21,0.23,0.24,0.25,0.27,0.28,0.29,0.31,0.32,0.33,0.35,0.36,0.37,0.39,0.40,0.41,0.43,0.44,0.45,0.47,0.48,0.49,0.51,0.52,0.53,0.55,0.56,0.57,0.59,0.60,0.61,0.63,0.64,0.65,0.67,0.68,0.69,0.71,0.72,0.73,0.75,0.76,0.77,0.79,0.80,0.81,0.83,0.84,0.85,0.87,0.88,0.89,0.91,0.92,0.93,0.95,0.96,0.97,0.99,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00};
        float temp_color_g[200]={0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.01,0.03,0.04,0.05,0.07,0.08,0.09,0.11,0.12,0.13,0.15,0.16,0.17,0.19,0.20,0.21,0.23,0.24,0.25,0.27,0.28,0.29,0.31,0.32,0.33,0.35,0.36,0.37,0.39,0.40,0.41,0.43,0.44,0.45,0.47,0.48,0.49,0.51,0.52,0.53,0.55,0.56,0.57,0.59,0.60,0.61,0.63,0.64,0.65,0.67,0.68,0.69,0.71,0.72,0.73,0.75,0.76,0.77,0.79,0.80,0.81,0.83,0.84,0.85,0.87,0.88,0.89,0.91,0.92,0.93,0.95,0.96,0.97,0.99,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00};
        float temp_color_b[200]={0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58,0.60,0.62,0.64,0.66,0.68,0.70,0.72,0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.00};

        for (int ix=0;ix<200;ix++) {
            color_r[ix]=temp_color_r[ix];
            color_g[ix]=temp_color_g[ix];
            color_b[ix]=temp_color_b[ix];
        }
    }
}
