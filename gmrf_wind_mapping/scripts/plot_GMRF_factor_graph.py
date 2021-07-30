#!/usr/bin/env python
"""
    Plots the GMRF factor-graph with the relations between the nodes
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image
import math


dataJ = "/home/jgmonroy/catkin_ws/src/olfaction/gas_distribution_mapping/gmrf_wind_mapping/scripts/gmrf_Jacobian.txt"
dataA = "/home/jgmonroy/catkin_ws/src/olfaction/gas_distribution_mapping/gmrf_wind_mapping/scripts/gmrf_Lambda.txt"
dataObs = "/home/jgmonroy/catkin_ws/src/olfaction/gas_distribution_mapping/gmrf_wind_mapping/scripts/gmrf_observations.txt"
occupancy_img = "/home/jgmonroy/catkin_ws/src/olfaction/gas_distribution_mapping/gmrf_wind_mapping/scripts/MAPIRlab_furniture_walls_wind.pgm"
occupancy_res = 0.01;

# Load occupancy gridmap (background), factor graph and observations
img = mpimg.imread(occupancy_img)
img = np.flipud(img)                #vertical flip

img2 = Image.open(occupancy_img)
width, height = img2.size
print("Occupancy_img is %u x %u px" %(width, height))
fileObs = open(dataObs,"r")
obs = fileObs.read().splitlines()
fileJ = open(dataJ,"r")
fileA = open(dataA,"r")
A = fileA.read().splitlines()
plt.ion()   #iterative mode plotting
f, (ax1, ax2) = plt.subplots(1, 2)


# Parse the file containing the Jacobian triplets
count = 0       #overal line counter in file
cells_x = 0     #num cells x
cells_y = 0     #num cells y
cell_size = 0   #m
Jcount = 0      #counter the elements for each row in the jacobian
L = []          #list of the elements of a Jacobian row (J_row Cell_idx J_value)


for line in fileJ:
    l = line.split()
    
    if count == 1:          # cells_x
        cells_x = int(l[3])
    elif count == 2:        # cells_y
        cells_y = int(l[3])
        N = cells_x * cells_y        
    elif count == 3:        # cell_size
        cell_size = float(l[3])
        print("GMRF with %u cells and %.2f (m) cell_size" % (N, cell_size) )        
        ax1.imshow(img,origin='lower',cmap='Greys_r')
        ax1.set_title("Occupancy map")
        #ax2.imshow(img,origin='lower',cmap='Greys_r', extent=[0,cells_x*cell_size,-0.05,cells_y*cell_size])
        ax2.imshow(img,origin='lower',cmap='Greys_r', extent=[0,width*occupancy_res,-0.05,height*occupancy_res])
    elif count > 3:
        
        if Jcount == 0:             # keep data
            L = l
            Jcount = Jcount +1
        elif int(l[0])==int(L[0]):  # append data (same J row) and continue
            L.extend(l)
            Jcount = Jcount +1
        else:                       # another row of the Jacobian! plot and restart counter
            #print(Jcount)            
            #print(L)
            if Jcount == 1:         # obstacle or observation factors
                if int(L[1])<N:
                    c1x = int(L[1]) % cells_x
                    c1y = int(math.floor(int(L[1])/cells_x))
                    if float(obs[int(L[0])]) != 0.0:
                        # Observation factor
                        print("new obs factor with value (%.2f,%.2f) at cell %lu" % (float(obs[int(L[0])]),float(obs[int(L[0])+1]),int(L[1]) ) )
                        ax2.plot([c1x*cell_size + cell_size/2],[c1y*cell_size + cell_size/2],'ro')
                    else:
                        #Obstacle factor (Wx=0)
                        ax2.plot([c1x*cell_size + cell_size/2],[c1y*cell_size + cell_size/2],'bo')
                else:
                    if float(obs[int(L[0])]) == 0.0:
                        #Obstacle factor (Wy=0)
                        c1x = (int(L[1])-N) % cells_x
                        c1y = int(math.floor((int(L[1])-N)/cells_x))
                        ax2.plot([c1x*cell_size+0.1 + cell_size/2],[c1y*cell_size + cell_size/2],'yo')
            elif Jcount == 2:       # regularization factor
                #Plot only for cell < N
                if int(L[1])<N:
                    c1x = int(L[1]) % cells_x
                    c1y = int(math.floor(int(L[1])/cells_x))
                    c2x = int(L[4]) % cells_x
                    c2y = int(math.floor(int(L[4])/cells_x))
                    ax2.plot([c1x*cell_size + cell_size/2, c2x*cell_size + cell_size/2],[c1y*cell_size + cell_size/2,c2y*cell_size + cell_size/2],'b')
                    #plt.tight_layout()
                    #plt.savefig('factor_graph.pdf')
                
            else:                   # mass-law factor
                print("Mass-law factor!")
                
            # Keep current triplet and start over
            L = l
            Jcount = 1
        #end-if
    count = count+1
    
    #end-if
#end for
ax2.set_title("GMRF factor graph")            
#plt.show()
print("factor graph created!")
file.close() 

"""
for line in file:    
    l = line.split()
    
    if count == 1:          # cells_x
        cells_x = int(l[3])
    elif count == 2:        # cells_y
        cells_y = int(l[3])
        N = cells_x * cells_y        
    elif count == 3:        # cell_size
        cell_size = float(l[3])
        print("GMRF with %u cells and %.2f (m) cell_size" % (N, cell_size) )        
        ax1.imshow(img,origin='lower',cmap='Greys_r')
        ax1.set_title("Occupancy map")
        ax2.imshow(img,origin='lower',cmap='Greys_r', extent=[0,cells_x*cell_size,-0.05,cells_y*cell_size])
    elif count > 2 and count < 10000:
        
        
        if int(l[1])<N: #only the Wx related factors
        
            # regularization factors
            if int(l[0])==idxJ and value==1 and float(l[2])==-1:
                #Plot            
                c1x = idx % cells_x
                c1y = int(math.floor(idx/cells_x))
                c2x = int(l[1]) % cells_x
                c2y = int(math.floor(int(l[1])/cells_x))
                #plot with cells
                #ax2.plot([c1x, c2x],[c1y,c2y],'bo')
                #ax2.plot([c1x, c2x],[c1y,c2y],'k')
                #plot with cell_size            
                #ax2.plot([c1x*cell_size + cell_size/2, c2x*cell_size + cell_size/2],[c1y*cell_size + cell_size/2,c2y*cell_size + cell_size/2],'bo')
                ax2.plot([c1x*cell_size + cell_size/2, c2x*cell_size + cell_size/2],[c1y*cell_size + cell_size/2,c2y*cell_size + cell_size/2],'b')
                
                #plt.tight_layout()
                #plt.savefig('factor_graph.pdf')
            
            # Null factor due to obstacles
            if int(l[0])!=idxJ and value==1:
                #Plot            
                c1x = idx % cells_x
                c1y = int(math.floor(idx/cells_x))                
                #ax2.plot([c1x*cell_size + cell_size/2, c2x*cell_size + cell_size/2],[c1y*cell_size + cell_size/2,c2y*cell_size + cell_size/2],'bo')
                ax2.plot([c1x*cell_size + cell_size/2],[c1y*cell_size + cell_size/2],'ro')        
        #update indexes
        idxJ = int(l[0])
        idx = int(l[1])
        value = float(l[2])
    count = count+1
#end-for    
ax2.set_title("GMRF factor graph")            
#plt.show()
print("factor graph created!")
file.close()
"""