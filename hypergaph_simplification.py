import skeleton
import matplotlib.pyplot as plt
import numpy as np

"""structure of an hyperedges should be like so : 
H = [(3, 2), (2, 14), (7,4)] 
tuple (i,j) : i corresponds to the vertex index, and j the index of said edge in edges_list. So the final list of pixels is edges_list[i][j]
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
            hyperedge_to_split = hyperedges[split_index].copy()
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
            merge_index = np.random.randint(0,len(hyperedges))
            edge1 = hyperedges[merge_index].copy()
            degree = degrees[merge_index]
            vertex = hyperedges[merge_index][0]

            # must find the hyperedges with a common vertex
            list = []
            for j in range(len(hyperedges)):
                if j != merge_index : 
                    hyperedge = hyperedges[j]
                    if hyperedge[0][0] == vertex or edges_list[hyperedge[-1][0]][hyperedge[-1][1]][-1] == vertex :
                        list.append(j)

            other_merge_index = np.random.choice(list)
            edge2 = hyperedges[other_merge_index].copy()

            if edge2[0] ==  vertex : #both edges start with the same vertex
                edge2.reverse()
            edge2 = edge2 + edge1
            del hyperedges[merge_index]
            del hyperedges[other_merge_index]
            del degrees[merge_index]
            del degrees[other_merge_index]
            hyperedges.append(edge2)
            degrees.append(degree)
    return hyperedges, degrees

