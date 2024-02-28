import numpy as np


##########################################################

#trianglur lattice of cells
#periodic driver: udot = b0 + u^2 until u = 1, then u = 0
#coupled periodic responders: vidot = bi + vi^2 + g0*(u-vi) + gSUM[(vj-vi)] 
# until vi = 1, then vi = 0, j subscript is nearest neighbors
#deviation score is calculated

##########################################################
#Parameters

row_size = 10
vmax = 1
vrs = 0
g0 = 0.05
b0 = 0.0025
T0 = 1/np.sqrt(b0)
g = 0.05

IC = 0

#probability that a node disconnected from all of its neighbors
p_node = 0.0

##############################################################

#create time array
t_start = 0
t_end = 8*T0
dt = 0.0001
ts = np.arange(t_start+dt,t_end,dt)

##############################################################

def QIF(row_size,vmax,vrs,g0,b0,g,ts,p_node):

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


    nodes = np.ones(num_of_cells)
    # turning off cells (nodes)
    for i in range(0,num_of_cells):
        rn = np.random.rand()
        if rn < p_node:
            nodes[i] = 0


    # creating arrays of links
    # note that extra links that are never used are added onto 
    # certain nodes. Doing this simplifies the indexing

    # left-right links
    r_lr = np.ones(num_of_cells)
    # bottom-right/top-left links
    r_brtl = np.ones((row_size - 1)*row_size)
    # bottom-left/top-right links
    r_bltr = np.ones((row_size - 1)*row_size)

    # relations between positions for links

    #l_l = r_lr[j-1]
    #l_r = r_lr[j]

    #For odd rows
    #l_otl = r_brtl[j-(row_size+1)]
    #l_otr = r_bltr[j-row_size]
    #l_obl = r_bltr[j]
    #l_obr = r_brtl[j]

    #For even rows
    #l_etl = r_brtl[j-row_size]
    #l_etr = r_bltr[j-(row_size-1)]
    #l_ebl = r_bltr[j]
    #l_ebr = r_brtl[j]

    # turning off links from off cells
    for j in range(0,num_of_cells):
        if nodes[j] == 0:
            #1st row
            #Note: First row is always odd and then they alternate
            if j  // row_size == 0:

                #1st column
                if j % row_size == 0:
                    r_lr[j] = 0
                    r_brtl[j] = 0

                #middle columns
                elif j % row_size != 0 and j % row_size != (row_size-1):
                    r_lr[j-1] = 0
                    r_lr[j] = 0
                    r_bltr[j] = 0
                    r_brtl[j] = 0

                #last column
                elif j % row_size == (row_size-1):
                    r_lr[j-1] = 0
                    r_bltr[j] = 0
                    r_brtl[j] = 0

            
            #middle rows
            elif j // row_size != 0 and j // row_size != (row_size - 1):
                
                #even rows
                if (j  // row_size) % 2 != 0:

                    #1st column
                    if j % row_size == 0:
                        r_brtl[j-row_size] = 0
                        r_bltr[j-(row_size-1)] = 0
                        r_lr[j] = 0

                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        r_lr[j] = 0
                        r_lr[j-1] = 0
                        r_brtl[j-row_size] = 0
                        r_bltr[j-(row_size-1)] = 0 
                        r_bltr[j] = 0
                        r_brtl[j] = 0

                    #last column
                    elif j % row_size == (row_size-1):
                        r_lr[j-1] = 0
                        r_brtl[j-row_size] = 0
                        r_bltr[j] = 0


                #odd rows
                else: #(j // row_size) % 2 == 0

                    #1st column
                    if j % row_size == 0:
                        r_bltr[j-row_size] = 0
                        r_bltr[j] = 0
                        r_lr[j] = 0

                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        r_lr[j] = 0
                        r_lr[j-1] = 0
                        r_brtl[j-(row_size+1)] = 0 
                        r_bltr[j-row_size] = 0
                        r_bltr[j] = 0
                        r_brtl[j] = 0
                        

                    #last column
                    elif j % row_size == (row_size-1):
                        r_lr[j-1] = 0
                        r_brtl[j-(row_size+1)] = 0 
                        r_bltr[j-row_size] = 0
                        r_bltr[j] = 0
                        r_brtl[j] = 0
                        

            #last row
            elif j // row_size == (row_size - 1):

                #even rows
                if (j // row_size) % 2 != 0:

                    #1st column
                    if j % row_size == 0:
                        r_brtl[j-row_size] = 0
                        r_bltr[j-(row_size-1)] = 0
                        r_lr[j] = 0

                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        r_lr[j-1] = 0
                        r_lr[j] = 0
                        r_brtl[j-row_size] = 0
                        r_bltr[j-(row_size-1)] = 0

                    #last column
                    elif j % row_size == (row_size-1):
                        r_brtl[j-row_size] = 0
                        r_lr[j-1] = 0

                #odd rows
                else: #(j // row_size) % 2 == 0

                    #1st column
                    if j % row_size == 0:
                        r_bltr[j-row_size] = 0
                        r_lr[j] = 0 

                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        r_lr[j-1] = 0
                        r_lr[j] = 0
                        r_brtl[j-(row_size+1)] = 0
                        r_bltr[j-row_size] = 0 

                    #last column
                    elif j % row_size == (row_size-1):
                        r_lr[j-1] = 0
                        r_brtl[j-(row_size+1)] = 0
                        r_bltr[j-row_size] = 0


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
            #l_l = r_lr[j-1]
            #v_r = vs[j+1]
            #l_r = r_lr[j]
            
            #For odd rows
            #v_otl = vs[j-(row_size+1)]
            #l_otl = r_brtl[j-(row_size+1)]
            #v_otr = vs[j-row_size]
            #l_otr = r_bltr[j-row_size]
            #v_obl = vs[j+(row_size-1)]
            #l_obl = r_bltr[j]
            #v_obr = vs[j+row_size]
            #l_obr = r_brtl[j]
            
            #For even rows
            #v_etl = vs[j-row_size]
            #l_etl = r_brtl[j-row_size]
            #v_etr = vs[j-(row_size-1)]
            #l_etr = r_bltr[j-(row_size-1)]
            #v_ebl = vs[j+row_size]
            #l_ebl = r_bltr[j]
            #v_ebr = vs[j+(row_size+1)]
            #l_ebr = r_brtl[j]

            #1st row
            #Note: First row is always odd and then they alternate
            if j  // row_size == 0:

                #1st column
                if j % row_size == 0:
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_r = vs[j+1]
                        l_r = r_lr[j]
                        v_obr = vs[j+row_size]
                        l_obr = r_brtl[j]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_obr*g*(v_obr[t] - v[t]))

                #middle columns
                elif j % row_size != 0 and j % row_size != (row_size-1):
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_l = vs[j-1]
                        l_l = r_lr[j-1]
                        v_r = vs[j+1]
                        l_r = r_lr[j]
                        v_obl = vs[j+(row_size-1)]
                        l_obl = r_bltr[j]
                        v_obr = vs[j+row_size]
                        l_obr = r_brtl[j]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_obl*g*(v_obl[t] - v[t]) + l_obr*g*(v_obr[t] - v[t]))

                #last column
                elif j % row_size == (row_size-1):
                    if v[t] > vmax:
                        v[t+1] = vrs
                    else:
                        v_l = vs[j-1]
                        l_l = r_lr[j-1]
                        v_obl = vs[j+(row_size-1)]
                        l_obl = r_bltr[j]
                        v_obr = vs[j+row_size]
                        l_obr = r_brtl[j]
                        v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_obl*g*(v_obl[t] - v[t]) + l_obr*g*(v_obr[t] - v[t]))


            #middle rows
            elif j // row_size != 0 and j // row_size != (row_size - 1):
                
                #even rows
                if (j  // row_size) % 2 != 0:

                    #1st column
                    if j % row_size == 0:
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_r = vs[j+1]
                            l_r = r_lr[j]
                            v_etl = vs[j-row_size]
                            l_etl = r_brtl[j-row_size]
                            v_etr = vs[j-(row_size-1)]
                            l_etr = r_bltr[j-(row_size-1)]
                            v_ebl = vs[j+row_size]
                            l_ebl = r_bltr[j]
                            v_ebr = vs[j+(row_size+1)]
                            l_ebr = r_brtl[j]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_etl*g*(v_etl[t] - v[t]) + l_etr*g*(v_etr[t] - v[t]) + l_ebl*g*(v_ebl[t] - v[t]) + l_ebr*g*(v_ebr[t] - v[t]))

                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_l = vs[j-1]
                            l_l = r_lr[j-1]
                            v_r = vs[j+1]
                            l_r = r_lr[j]
                            v_etl = vs[j-row_size]
                            l_etl = r_brtl[j-row_size]
                            v_etr = vs[j-(row_size-1)]
                            l_etr = r_bltr[j-(row_size-1)]
                            v_ebl = vs[j+row_size]
                            l_ebl = r_bltr[j]
                            v_ebr = vs[j+(row_size+1)]
                            l_ebr = r_brtl[j]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_etl*g*(v_etl[t] - v[t]) + l_etr*g*(v_etr[t] - v[t]) + l_ebl*g*(v_ebl[t] - v[t]) + l_ebr*g*(v_ebr[t] - v[t]))

                    #last column
                    elif j % row_size == (row_size-1):
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_l = vs[j-1]
                            l_l = r_lr[j-1]
                            v_etl = vs[j-row_size]
                            l_etl = r_brtl[j-row_size]
                            v_ebl = vs[j+row_size]
                            l_ebl = r_bltr[j]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_etl*g*(v_etl[t] - v[t]) + l_ebl*g*(v_ebl[t] - v[t]))


                #odd rows
                else: #(j // row_size) % 2 == 0

                    #1st column
                    if j % row_size == 0:
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_r = vs[j-1]
                            l_r = r_lr[j-1]
                            v_otr = vs[j-row_size]
                            l_otr = r_bltr[j-row_size]
                            v_obr = vs[j+row_size]
                            l_obr = r_brtl[j]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_otr*g*(v_otr[t] - v[t]) + l_obr*g*(v_obr[t] - v[t]))

                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_l = vs[j-1]
                            l_l = r_lr[j-1]
                            v_r = vs[j-1]
                            l_r = r_lr[j-1]
                            v_otl = vs[j-(row_size+1)]
                            l_otl = r_brtl[j-(row_size+1)]
                            v_otr = vs[j-row_size]
                            l_otr = r_bltr[j-row_size]
                            v_obl = vs[j+(row_size-1)]
                            l_obl = r_bltr[j]
                            v_obr = vs[j+row_size]
                            l_obr = r_brtl[j]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_otl*g*(v_otl[t] - v[t]) + l_otr*g*(v_otr[t] - v[t]) + l_obl*g*(v_obl[t] - v[t]) + l_obr*g*(v_obr[t] - v[t]))

                    #last column
                    elif j % row_size == (row_size-1):
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_l = vs[j-1]
                            l_l = r_lr[j-1]
                            v_otl = vs[j-(row_size+1)]
                            l_otl = r_brtl[j-(row_size+1)]
                            v_otr = vs[j-row_size]
                            l_otr = r_bltr[j-row_size]
                            v_obl = vs[j+(row_size-1)]
                            l_obl = r_bltr[j]
                            v_obr = vs[j+row_size]
                            l_obr = r_brtl[j]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_otl*g*(v_otl[t] - v[t]) + l_otr*g*(v_otr[t] - v[t]) + l_obl*g*(v_obl[t] - v[t]) + l_obr*g*(v_obr[t] - v[t]))


            #last row
            elif j // row_size == (row_size - 1):

                #even rows
                if (j // row_size) % 2 != 0:

                    #1st column
                    if j % row_size == 0:
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_r = vs[j+1]
                            l_r = r_lr[j]
                            v_etl = vs[j-row_size]
                            l_etl = r_brtl[j-row_size]
                            v_etr = vs[j-(row_size-1)]
                            l_etr = r_bltr[j-(row_size-1)]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_etl*g*(v_etl[t] - v[t]) + l_etr*g*(v_etr[t] - v[t]))


                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_l = vs[j-1]
                            l_l = r_lr[j-1]
                            v_r = vs[j+1]
                            l_r = r_lr[j]
                            v_etl = vs[j-row_size]
                            l_etl = r_brtl[j-row_size]
                            v_etr = vs[j-(row_size-1)]
                            l_etr = r_bltr[j-(row_size-1)]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_etl*g*(v_etl[t] - v[t]) + l_etr*g*(v_etr[t] - v[t]))


                    #last column
                    elif j % row_size == (row_size-1):
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_l = vs[j-1]
                            l_l = r_lr[j-1]
                            v_etl = vs[j-row_size]
                            l_etl = r_brtl[j-row_size]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_etl*g*(v_etl[t] - v[t]))


                #odd rows
                else: #(j // row_size) % 2 == 0

                    #1st column
                    if j % row_size == 0:
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_r = vs[j-1]
                            l_r = r_lr[j-1]
                            v_otr = vs[j-row_size]
                            l_otr = r_bltr[j-row_size]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_otr*g*(v_otr[t] - v[t]))


                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_l = vs[j-1]
                            l_l = r_lr[j-1]
                            v_r = vs[j-1]
                            l_r = r_lr[j-1]
                            v_otl = vs[j-(row_size+1)]
                            l_otl = r_brtl[j-(row_size+1)]
                            v_otr = vs[j-row_size]
                            l_otr = r_bltr[j-row_size]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_r*g*(v_r[t] - v[t]) + l_otl*g*(v_otl[t] - v[t]) + l_otr*g*(v_otr[t] - v[t]))


                    #last column
                    elif j % row_size == (row_size-1):
                        if v[t] > vmax:
                            v[t+1] = vrs
                        else:
                            v_l = vs[j-1]
                            l_l = r_lr[j-1]
                            v_otl = vs[j-(row_size+1)]
                            l_otl = r_brtl[j-(row_size+1)]
                            v_otr = vs[j-row_size]
                            l_otr = r_bltr[j-row_size]
                            v[t+1] = v[t] + dt*(b + v[t]**2 + g0*(u[t] - v[t]) + l_l*g*(v_l[t] - v[t]) + l_otl*g*(v_otl[t] - v[t]) + l_otr*g*(v_otr[t] - v[t]))


      
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
num_of_trials = 15

DSs = np.array([])

for i in range(0,num_of_trials):
    vs = QIF(row_size,vmax,vrs,g0,b0,g,ts,p_node)

    DS = Deviation_Score(vs, T0, dt, t_end)
    DSs = np.append(DSs, DS)


