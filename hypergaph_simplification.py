import skeleton
import matplotlib.pyplot as plt
import numpy as np

"""structure of an hyperedges should be like so : 
H = [(3, 2), (2, 14), (14,7)] 
tuple (i,j) : where i is the start point and j the end point. (i,j) can be found in edges_list[i] such that edge[-1] == j
"""

def fidelity_energy(skeleton, hyperedges : list) :
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
    width = len(skeleton[0])
    energy = 0
    for hyperedge in hyperedges : 
        i, j = hyperedge
        edge = 0
        edge = hyperedges[i][j]
        P0 = edge[0]
        P3 = edge[-1]
        control_points = skeleton.get_control_points_bezier(edge, P0, P3, width)
        energy += skeleton.fitting_error(edge, control_points, width)
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
    return energy

def trade_off_energy(skeleton, hyperedges : list, mu : float, alpha : float, degrees : list) :
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
    return (1-alpha)*fidelity_energy(skeleton, hyperedges) + alpha*simplificity_energy(degrees, mu)


def merge_and_split(hyperedges: list, degrees : list, edges_list, max_changes = 15) :
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

    if mode == 0 : #split
        #must choose nb_changes random hyperedges and split them at a random place.
        for i in range(nb_changes):
            split_index = np.random.randint(0,len(hyperedges))
            #must get the values to use them further on
            degree = degrees[split_index].copy()
            hyperedge_to_split = hyperedges[split_index].copy() #list of tuples
            if len(hyperedge_to_split)>1:
                #split the edge
                split_point = np.random.randint(1, len(hyperedge_to_split)-1)
                new_edge1 = hyperedge_to_split[0:split_point]
                new_edge2 = hyperedge_to_split[split_point:]

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

            vertex = hyperedges[merge_index][0] #merging vertex
            # must find the hyperedges starting or ending by vertex
            list = []
            for j in range(len(hyperedges)):
                if j != merge_index : 
                    hyperedge = hyperedges[j]
                    if hyperedge[0][0] == vertex or hyperedge[-1][-1] == vertex :
                        list.append(j)

            other_merge_index = np.random.choice(list)
            while hyperedges[other_merge_index][0] == hyperedge1[0] or hyperedges[other_merge_index][-1][0] == hyperedge1[0][1] :
                #for convenience, I'd rather not handle the case where the hyperedges are included one into the other one
                other_merge_index = np.random.choice(list)

            hyperedge2 = hyperedges[other_merge_index].copy()

            if hyperedge2[0][0] ==  vertex :
                hyperedge2.reverse() #both edges start with the same vertex
                for e in range(len(hyperedge2)) :
                    i, j = hyperedge2[e]
                    hyperedge2[e] = (j,i)


            hyperedge2 = hyperedge2 + hyperedge1
            del hyperedges[merge_index]
            del hyperedges[other_merge_index]
            del degrees[merge_index]
            del degrees[other_merge_index]
            hyperedges.append(hyperedge2)
            degrees.append(degree)
    return hyperedges, degrees


def degree_switch(degrees, proba):
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

def convert_to_hypergraph(edges_list):
    hypergraph = []
    for edge_list in edges_list : 
        for edge in edge_list : 
            i = edge[0]
            j = edge[-1]
            hypergraph.append([(i,j)])
    return hypergraph

def modified_metropolis_hastings(hyperedges, edges_list, degrees, junctions_indices, Tinit : float = 1, flip_proba = 0.3):
    """
    
    Parameters : 
    --------------
    
    Returns : 
    --------------
    
    """
    if len(junctions_indices) < 1 : 
        raise Exception("graph is empty")
        exit(0)
    hyperedges = convert_to_hypergraph(edges_list)
    T = Tinit
    V = len(junctions_indices)
    C = 0.999**(1/V)
    T_end = C**10000
    while T > T_end : 
        perturbation_operator = np.random.randint(0, 3)
        if perturbation_operator == 0 : 
            hyperedges, degrees = merge_and_split(hyperedges, degrees, edges_list)
        elif perturbation_operator == 1 : 
            degrees = degree_switch(degrees, flip_proba)
            