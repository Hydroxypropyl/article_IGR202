
import skeleton
import matplotlib.pyplot as plt
import numpy as np

"""structure of an hyperedges should be like so : 
H = [(3, 2), (2, 14), (14,7)] 
tuple (i,j) : where i is the start point and j the end point. (i,j) can be found in edges_list[i] such that edge[-1] == j
"""

def fidelity_energy(width, hyperedges : list, junctions_indices, edges_list, degrees : list) :
    """
    Takes as input 
    
    Parameters :
    -------------
    `skeleton` : 
    `hyperedges` : the hyperedges. Coordinates are in flattened image.
    
    Returns : 
    -------------
    The fidelity energy as the sum of the fitting energy on all the hyperedges.
    
    """ 
    energy = 0
    for h in range(len(hyperedges)) : 
        degree = degrees[h]
        hyperedge = hyperedges[h]
        edge = [hyperedge[0][0]]
        for edge_extremities in hyperedge :
            P0 = edge_extremities[0]
            P3 = edge_extremities[-1]
            i = junctions_indices.index(P0)
            for e in edges_list[i]:
                if e[-1] == P3 :
                    edge = edge + e[1:]
                    break
            if edge == 0 : 
                raise Exception("issue while finding the corresponding edge ({},{}).".format(P0,P3))
                exit(0)
        control_points = skeleton.get_control_points_bezier(edge, edge[0], edge[-1], width, degree)
        energy += skeleton.fitting_error(edge, control_points, width)
    #print("fidelity energy :", energy)
    return energy


def simplificity_energy(degrees : list, mu : float) :
    """
    
    Parameters : 
    ---------------
    
    `deg` : degree of the bezier curve. Must be 1, 2 or 3.
    `mu` : model parameter. Higher value increases the presence of straight lines.
    
    Returns : 
    ---------------

    """
    energy = 0
    for degree in degrees :
        energy += 1 + mu*degree
    #print("simplicity energy :", energy)
    return energy

def trade_off_energy(width, junctions_indices, hyperedges : list, mu : float, alpha : float, degrees : list, edges_list) :
    """"
    Computes the trade-off energy using the simplicity and fidelity energy.
    
    Parameters :
    ----------------
    `skeleton`:
    `hyperedges`:
    `mu`:
    `alpha`:
    `degrees` : degrees of the bezier curves for the corresponding hyperedges.
    
    Returns : 
    ----------------
    
    """
    return (1-alpha)*fidelity_energy(width, hyperedges, junctions_indices, edges_list, degrees) + alpha*simplificity_energy(degrees, mu)


def merge_and_split(hyperedges: list, degrees : list, edges_list : list, max_changes = 15) -> tuple[list, list] :
    """
    Randomly merges or splits a random number of hyperedges.
    
    Parameters:
    --------------
    `hyperedges`: list of hyperedges.
    `degrees`:
    `edges_list`:
    `max_changes` : upper bound for the random generation of the number of changes to do.
    
    Returns :
    --------------
    `hyperedges`: the updated hyperedges
    `degrees` : updated degrees
    """

    mode = np.random.randint(0,2)
    nb_changes = np.random.randint(0,max_changes)
    #print("mode : ", mode)
    if mode == 0 : #split
        #must choose nb_changes random hyperedges and split them at a random place.
        for i in range(nb_changes):
            split_index = np.random.randint(0,len(hyperedges))
            #must get the values to use them further on
            degree = degrees[split_index]
            hyperedge_to_split = hyperedges[split_index].copy() #list of tuples
            if len(hyperedge_to_split)>1:
                #split the edge
                split_point = np.random.randint(1, len(hyperedge_to_split))
                new_edge1 = hyperedge_to_split[0:split_point]
                new_edge2 = hyperedge_to_split[split_point:]
                #print("split point {}, new edges : {} and {}".format(split_point, new_edge1, new_edge2))

                del hyperedges[split_index]
                del degrees[split_index]
                #update the hyperedges and the degrees
                hyperedges.append(new_edge1)
                degrees.append(degree)
                hyperedges.append(new_edge2)
                degrees.append(degree)                     
    else : # merge
        for i in range(nb_changes):
            #choose the hyperedge to merge
            merge_index = np.random.randint(0,len(hyperedges))
            hyperedge1 = hyperedges[merge_index].copy()

            degree = degrees[merge_index]

            vertex = hyperedges[merge_index][0][0] #merging vertex
            # must find the hyperedges starting or ending by vertex
            #print("vertex is {}, merge index is {}, hyperedge1 {}".format(vertex, merge_index, hyperedge1))
            list = []
            for j in range(len(hyperedges)):
                if j != merge_index : 
                    hyperedge = hyperedges[j]
                    #print("candidate : {}".format(hyperedge))
                    if hyperedge[0][0] == vertex or hyperedge[-1][-1] == vertex :
                        #print("found a candidate !")
                        list.append(j)
            if len(list) > 0 : 
                other_merge_index = np.random.choice(list)
                while hyperedges[other_merge_index][0] == hyperedge1[0] or hyperedges[other_merge_index][-1][0] == hyperedge1[0][1] :
                    list.remove(other_merge_index)
                    #for convenience, I'd rather not handle the case where the hyperedges are included one into the other one
                    if len(list) == 0 :
                        break
                    other_merge_index = np.random.choice(list)
                if len(list)> 0 : 
                    hyperedge2 = hyperedges[other_merge_index].copy()
                    #print("I've chosen candidate ", hyperedge2)

                    if hyperedge2[0][0] ==  vertex :
                        hyperedge2.reverse() #both edges start with the same vertex
                        for e in range(len(hyperedge2)) :
                            i, j = hyperedge2[e]
                            hyperedge2[e] = (j,i)


                    hyperedge2 = hyperedge2 + hyperedge1
                    if merge_index > other_merge_index :
                        del hyperedges[merge_index]
                        del hyperedges[other_merge_index]
                        del degrees[merge_index]
                        del degrees[other_merge_index]
                    else : 
                        del hyperedges[other_merge_index]
                        del hyperedges[merge_index]
                        del degrees[other_merge_index]
                        del degrees[merge_index]

                    hyperedges.append(hyperedge2)
                    degrees.append(degree)
    #print("final hyperedges {} \n final degrees{}".format(hyperedges, degrees))
    return hyperedges, degrees


def degree_switch(degrees : list, proba : float) -> list:
    """
    For each degree in `degrees`, there is a probibility inferior to `proba` that its value will be recomputed, by being randomly generated in {1, 2, 3}.
    
    Parameters : 
    --------------
    `degrees` : the list of degrees for all the bezier curves associated to the hyperedges.
    `proba` : the flip probability. Must be in [0, 1]

    Returns :
    --------------
    `degrees`

    """
    for d in range(len(degrees)) : 
        if np.abs(proba-1) < 10e-7 :
            degrees[d] = np.random.randint(1, 4)
        else : 
            flip = np.random.random(1)
            if flip <= proba : 
                degrees[d] = np.random.randint(1, 4)
    return degrees
                

def overlap_and_dissociation(hyperedges : list, max_changes = 15) :
    """
    
    Parameters : 
    ---------------
    `hyperedges` : 
    `max_changes` : 
    
    Returns :
    ---------------
    """
    mode = np.random.randint(0,2)
    nb_changes = np.random.randint(0,max_changes)
    if mode == 0 : #overlap
        #must choose nb_changes random hyperedges and split them at a random place.
        for i in range(nb_changes):
            pass
    else : # dissociate
        pass

def convert_to_hyperedges(edges_list : list, default_degree = 3) -> tuple[list, list]:
    """
    
    Parameters : 
    --------------
    
    Returns : 
    --------------
    `hyperedges`:
    `degrees`:
    """
    hyperedges = []
    for edge_list in edges_list : 
        for edge in edge_list : 
            i = edge[0]
            j = edge[-1]
            hyperedges.append([(i,j)])
    for hyperedge in hyperedges : 
        i,j = hyperedge[0]
        if [(j,i)] in hyperedges :
            hyperedges.remove([(j,i)])
    degrees = [default_degree for k in range(len(hyperedges))]
    return hyperedges, degrees


def hyperedges_fit_bezier(edges_list : list, hyperedges : list, degrees : list, width, junctions_indices):
    x = []
    y = []
    for i in range(len(hyperedges)) :
        degree = degrees[i]
        hyperedge = hyperedges[i]
        edge = [hyperedge[0][0]]
        for pair in hyperedge :
            start = pair[0]
            end = pair[1]
            index = junctions_indices.index(start)
            for e in edges_list[index]:
                if e[-1] == end :
                    edge = edge + e[1:]
                    break
            if edge == 0 : 
                raise Exception("issue while finding the corresponding edge ({},{}).".format(start,end))
                exit(0)
        control_points = skeleton.get_control_points_bezier(edge, edge[0], edge[-1], width, degree)
        norm = len(edge)
        for i in range(norm) :
            tp =i/norm
            bezier = skeleton.bezier_curve(tp, control_points)
            x.append(bezier[0])
            y.append(bezier[1])
    return x, y


def modified_metropolis_hastings(skeleton, edges_list : list, junctions_indices : list, Tinit : float = 1, flip_proba = 0.3, mu = 0.4, alpha = 0.05):
    """
    
    Parameters : 
    --------------
    
    Returns : 
    --------------
    `hyperedges` : 
    `degrees` :
    
    """
    if len(junctions_indices) < 1 : 
        raise Exception("graph is empty")
        exit(0)
    hyperedges, degrees = convert_to_hyperedges(edges_list)
    T = Tinit
    V = len(junctions_indices)
    C = 0.999**(1/V)
    T_end = C**25000
    width = len(skeleton[0])
    i = 1
    while T > T_end : 
        #print(hyperedges)
        energy = trade_off_energy(width, junctions_indices, hyperedges, mu, alpha, degrees, edges_list)
        perturbation_operator = np.random.randint(0,2)#np.random.randint(0, 3)
        #print("perturbation operator : {}".format(perturbation_operator))
        potential_hyperedges = 0
        potential_degrees = 0
        if perturbation_operator == 0 : 
            #print("im in mode 0")
            potential_hyperedges, potential_degrees = merge_and_split(hyperedges.copy(), degrees.copy(), edges_list)
        elif perturbation_operator == 1 : 
            #print("I'm in mode 1")
            potential_degrees = degree_switch(degrees.copy(), flip_proba)
            potential_hyperedges = hyperedges.copy()
            #print("debug")
        new_energy = trade_off_energy(width,junctions_indices, potential_hyperedges, mu, alpha, potential_degrees, edges_list)
        p = np.random.random(1)
        #print("T = {} and energies = {}, new {} ".format(T, energy, new_energy))
        if (energy-new_energy)/T >-30000 : 
            if p < np.exp((energy-new_energy)/T) : 
                #print("update")
                hyperedges, degrees = potential_hyperedges, potential_degrees
        T = T*C
        if i//1000*1000 == i :
            print("{}th iteration. There are {} hyperedges. Energy is {}.\n".format(i, len(hyperedges), energy))
        i += 1
    return hyperedges, degrees

test = True
if test : 
    name = 'boule'
    img, edges_list, junctions_indices = skeleton.test(name)
    width = len(img[0])
    print("computing hyperedges...")
    hyperedges, degrees = modified_metropolis_hastings(img, edges_list, junctions_indices)
    print("done.\n")
    print("number of edges : {}; number of hyperedges : {}".format(len(edges_list), len(hyperedges)))
    print("computing bezier...")
    x, y = hyperedges_fit_bezier(edges_list, hyperedges, degrees, width, junctions_indices)
    print("done.\n")
    fig = plt.figure(figsize=(15,15))
    plt.plot(y, x, color='black', marker = 'o', linestyle ='none', ms = 1)
    ax = plt.gca()
    ax.set_xlim(0,len(img[0]))
    ax.set_ylim(len(img),0)
    ax.set_aspect('equal', adjustable='box')
    plt.title("final bezier curves : {} hyperedges".format(len(hyperedges)))
    plt.savefig("results/"+name+"_final_result.png")
    plt.show()

debugging_merge_split = False
if debugging_merge_split : 
    hyperedges, degrees = merge_and_split([[(2, 3), (5, 1)], [(2, 3), (2, 2), (3, 1)], [(2,1)]], [2, 2, 2], [0])
    print("hyperedges : {} \n degrees : {}".format(hyperedges, degrees))
