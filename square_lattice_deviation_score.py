import numpy as np

##########################################################

#square lattice of cells
#periodic driver: udot = b0 + u^2 until u = 1, then u = 0
#coupled periodic responders: vidot = bi + vi^2 + g0*(u-vi) + gSUM[(vj-vi)] 
# until vi = 1, then vi = 0, j subscript is nearest neighbors
#deviation score is calculated

##########################################################
#Parameters

row_size = 18
vmax = 1
vrs = 0
g0 = 0.05
b0 = 0.0025
T0 = 1/np.sqrt(b0)
g = 0.05

IC = 0

#probability that an edge (connection) between nodes exists
p_couple = 0.0

##############################################################

#create time array
t_start = 0
t_end = 8*T0
dt = 0.0001
ts = np.arange(t_start+dt,t_end,dt)

##############################################################

def QIF(row_size,vmax,vrs,g0,b0,g,ts,p_couple):

    num_of_cells = row_size**2
    #creating an array of all vs, besides the initial condition, 0s are added as place holders
    vs = np.zeros((num_of_cells,np.size(ts)))
    for i in range(0,num_of_cells):
        vs[i,0] = IC


    #creating array of b values
    bs = np.random.normal(b0,0.1*b0,num_of_cells)
    for i in range(0,len(bs)):
        if bs[i] <= 0:
            bs[i] = 0.1*b0


    #Creating arrays to see if links are there or not
    #left-right links: note that an extra link that is never used is added onto 
    # the nodes in the last column. Doing this simplifies the indexing
    num_lr_links = num_of_cells
    r_lr = np.random.rand(num_lr_links)
    for ri in range(np.size(r_lr)):
        if r_lr[ri] < p_couple:
            r_lr[ri] = 1
        else:
            r_lr[ri] = 0

    #up-down links
    num_ud_links = row_size*(row_size-1)
    r_ud = np.random.rand(num_ud_links)
    for ri in range(np.size(r_ud)):
        if r_ud[ri] < p_couple:
            r_ud[ri] = 1
        else:
            r_ud[ri] = 0


    #periodic driver:
    u = np.array([IC])

    #generate dynamics
    for t in range(0,np.size(ts)-1):

        if u[t] > vmax:
            u = np.append(u,vrs)
        else:
            u_pp = u[t] + dt*(b0 + (u[t])**2)
            u = np.append(u,u_pp)

        for j in range(0,num_of_cells):

            v = vs[j]
            b = bs[j]

            #v_l = vs[j-1]
            #v_r = vs[j+1]
            #v_u = vs[j-row_size]
            #v_d = vs[j+row_size]

            #l_l = r_lr[j-1]
            #l_r = r_lr[j]
            #l_u = r_ud[j-row_size]
            #l_d = r_ud[j]

            #1st row
            if j // row_size == 0:

                #1st column
                if j % row_size == 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_r = vs[j+1]
                        v_d = vs[j+row_size]
                        l_r = r_lr[j]
                        l_d = r_ud[j]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_d*g*(v_d[t] - v[t]))

                    
                #middle columns
                elif j % row_size != 0 and (j+1) % row_size != 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_l = vs[j-1]
                        v_r = vs[j+1]
                        v_d = vs[j+row_size]
                        l_l = r_lr[j-1]
                        l_r = r_lr[j]
                        l_d = r_ud[j]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_d*g*(v_d[t] - v[t]))
                    
                    
                #last column
                elif (j+1) % row_size == 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_l = vs[j-1]
                        v_d = vs[j+row_size]
                        l_l = r_lr[j-1]
                        l_d = r_ud[j]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_d*g*(v_d[t] - v[t]))
                    
                
                
            #middle rows
            elif j // row_size != 0 and j // row_size != (row_size - 1):
                
                #1st column
                if j % row_size == 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_r = vs[j+1]
                        v_u = vs[j-row_size]
                        v_d = vs[j+row_size]
                        l_r = r_lr[j]
                        l_u = r_ud[j-row_size]
                        l_d = r_ud[j]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_u*g*(v_u[t] - v[t]) + l_d*g*(v_d[t] - v[t]))
                       
                        
                #middle columns
                elif j % row_size != 0 and (j+1) % row_size != 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_l = vs[j-1]
                        v_r = vs[j+1]
                        v_u = vs[j-row_size]
                        v_d = vs[j+row_size]
                        l_l = r_lr[j-1]
                        l_r = r_lr[j]
                        l_u = r_ud[j-row_size]
                        l_d = r_ud[j]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_u*g*(v_u[t] - v[t]) + l_d*g*(v_d[t] - v[t]))
                    
                    
                #last column
                elif (j+1) % row_size == 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_l = vs[j-1]
                        v_u = vs[j-row_size]
                        v_d = vs[j+row_size]
                        l_l = r_lr[j-1]
                        l_u = r_ud[j-row_size]
                        l_d = r_ud[j]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_u*g*(v_u[t] - v[t]) + l_d*g*(v_d[t] - v[t]))
                    
                        
                        
                    
            #last row
            elif j // row_size == (row_size - 1):
                
                #1st column
                if j % row_size == 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_r = vs[j+1]
                        v_u = vs[j-row_size]
                        l_r = r_lr[j]
                        l_u = r_ud[j-row_size]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_u*g*(v_u[t] - v[t]))
                    
                    
                #middle columns
                elif j % row_size != 0 and (j+1) % row_size != 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_l = vs[j-1]
                        v_r = vs[j+1]
                        v_u = vs[j-row_size]
                        l_l = r_lr[j-1]
                        l_r = r_lr[j]
                        l_u = r_ud[j-row_size]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_u*g*(v_u[t] - v[t]))

                    
                    
                #last column
                elif (j+1) % row_size == 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_l = vs[j-1]
                        v_u = vs[j-row_size]
                        l_l = r_lr[j-1]
                        l_u = r_ud[j-row_size]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_u*g*(v_u[t] - v[t]))
                      
        
    return vs


############################################################################################
def Deviation_Score(vs,T0,dt,t_end):
    vbar = np.array([])

    for j in range(0,np.size(vs[1])):
        v = 0 
        for i in range(0,np.size(vs,axis=0)):
            v += vs[i,j]
        vbar = np.append(vbar,(v/np.size(vs,axis=0)))
    
    diffs = np.array([])
    for i in range(0,np.size(vs,axis=0)):
        diffs = np.append(diffs,(dt/T0)*(T0/t_end)*(sum(abs(vs[i]-vbar)))) 

    deviation_score = np.average(diffs)
    
    return deviation_score

########################################################################################
#Running simulations
num_of_trials = 30

DSs = np.array([])

for i in range(0,num_of_trials):
    vs = QIF(row_size,vmax,vrs,g0,b0,g,ts,p_couple)

    DS = Deviation_Score(vs, T0, dt, t_end)
    DSs = np.append(DSs, DS)



