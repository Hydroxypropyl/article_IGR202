import matplotlib.pyplot as plt
from matplotlib.image import imsave
import numpy as np
from scipy import ndimage
import cv2 as cv
import networkx as nx

from skimage.morphology import dilation
from skimage.morphology import rectangle
from skimage.morphology import diamond
import skimage.io

def get_skeleton(img) :
    """
    Parameters : 
    -----------------
    `img` : the image for the skeleton extraction. Must be black lines on white background.
    
    Returns : 
    -----------------
    `skeleton` : same shape as `img`, with values in {0, 255}. Pixels belonging to the skeleton are white.
    """
    classes = list(set(img.flatten().tolist()))
    print("There are {} classes".format(len(classes)))
    masks = []
    print("Computing masks...")
    for nr_class in classes :
        mask = 1*(img==nr_class)
        masks.append(mask)
    print("done.")
    superposition = np.zeros(img.shape)
    print("Computing intersections and skeleton")
    for nr_class in range(len(classes)):
        print("executing for class {}. {} elements to consider".format(nr_class, len(classes)-nr_class))
        footprint = rectangle(3,3)
        img1 = dilation(masks[nr_class], footprint)
        for k in range(nr_class+1, len(classes)) :
            img2 =masks[k] 
            inter = 1*((img1 + img2) == 2)

            superposition += inter
    return 255*(superposition > 0)


def open_curves(skeleton, img):
    pass
    ##TODO

def is_junction(skeleton, pixel : tuple) :
    """
    Checks whether a certain pixel is a junction node.
    
    Parameters :
    --------------
    `skeleton` : image in black and white, [0, 255]
    `pixel` : tuple containing the x and y coordinates of the pixel
    
    Returns : 
    --------------
    `True` if `pixel` is a T, X or Y junction, `False` otherwise."""
    x = pixel[0]
    y = pixel[1]
    color = skeleton[x,y]
    
    #check 4-connexity : if 3 neighbours or more : junctions
    if skeleton[x-1,y] + skeleton[x, y-1] + skeleton[x,y+1] + skeleton[x+1, y] >= 3*color :
        return True
    else : 
        return False


def get_junctions(skeleton) :
    """
    
    Parameters : 
    --------------
    `skeleton` : image in black and white, [0, 255]

    Returns :
    --------------
    `junction_coordinates` : the coordinates of the junctions in the image.
    """
    junctions_indices = []
    height = len(skeleton)
    width = len(skeleton[0])
    for i in range(height):
        for j in range(width):
            if skeleton[i,j] > 0 : 
                bool = is_junction(skeleton, (i,j))
                if bool == True :
                    #print("junction at {},{}".format(i,j))
                    junctions_indices.append((i,j))
    return junctions_indices


def neighbours_8_connexity(img, p_i : int, p_j : int):
    """
    Computes the neighbours in 8-connexity in the image. 
    Img must be black and white and neighbours are those of same color as pixel p.
    Parameters:
    -----------
    `img` : array of integers
    `p_i`, `p_j` : integer coordinates of pixel of interest

    Returns : 
    -----------
    A list of tuples containing the coordinates of the neighbours in `img`

    """
    neighbours = []
    #print(img[p_i-1:p_i+2, p_j-1:p_j+2])
    candidates = [(p_i-1, p_j), (p_i, p_j-1), (p_i, p_j+1), (p_i+1, p_j), (p_i-1, p_j-1), (p_i-1, p_j+1), (p_i+1, p_j-1), (p_i+1, p_j+1)]
    for pixel in candidates : 
        if img[pixel[0], pixel[1]] == 255 :
                    neighbours.append((pixel[0],pixel[1]))
    #if len(neighbours) == 0 : 
     #   candidates = [(p_i-1, p_j-1), (p_i-1, p_j+1), (p_i+1, p_j-1), (p_i+1, p_j+1)]
      #  for pixel in candidates : 
       #     if img[pixel[0], pixel[1]] == 255 :
        #        neighbours.append((pixel[0],pixel[1]))

    if len(neighbours) == 0 :
        raise Exception("pixels has no neighbours???")
    return neighbours




def get_neighbours_and_edges(skeleton, junction_coordinates : list) -> tuple[list, list]:
    """
    Computes the edges determined by white pixels in `skeleton` between junctions.
    
    Parameters : 
    -------------------
    `skeleton` : image in black and white, [0, 255]
    `junction_coordinates` : list containing tuples of junctions coordinates in `skeleton`.
    
    Returns : 
    -------------------
    `junction_neighbours` : same length as `junction_coordinates`. 
    ith element contains the list of junctions connected by an edge to junction[i].
    Neighbours correspond to their index in `junction_coordinates`.
    `edges_list` : a list of same length as `junction_coordinates`, where edges_list[i] is the list of edges(all the pixels) for junction_coordinates[i]. 
    Edges contain the coordinates of the pixels.
    """

    edges_list =[]
    junction_neighbours=[]
    i = 1
    for junction in junction_coordinates : 
        print("junction {} out of {}".format(i, len(junction_coordinates)))
        i += 1
        #get neighbours : 
        edges = []
        junction_neighbours = neighbours_8_connexity(skeleton, junction[0], junction[1])
        nb_potential_edges = len(junction_neighbours)
        j = 1
        for n in range(nb_potential_edges) : 
            visited = np.zeros(skeleton.shape)
            edge = [junction, junction_neighbours[n]]
            print("edge {} out of {}".format(j, len(junction_neighbours)))
            j += 1
            not_edge = False
            last_neighbour = edge[-2]
            while is_junction(skeleton, edge[-1]) == False and not_edge == False :
                #find neighbour of latest pixel of edge
                
                pixel = edge[-1]
                #print("edge length : ", len(edge))
                neighbours = neighbours_8_connexity(skeleton, pixel[0], pixel[1]) # one of them should be edge[-2]
                if len(neighbours) < 2 : 
                    #must be an artifact
                    print("{} is in edge but has {} neighbours".format(pixel, len(neighbours)))
                    print(skeleton[pixel[0]-1:pixel[0]+2, pixel[1]-1:pixel[1]+2])
                    not_edge = True
                    
                else : 
                    print(neighbours, pixel, "last", last_neighbour)
                    print(skeleton[pixel[0]-1:pixel[0]+2, pixel[1]-1:pixel[1]+2])
                    for nbr in neighbours :
                        if nbr != last_neighbour and visited[nbr[0], nbr[1]] == 0:
                            #print(nbr, edge[-2])
                            edge.append(nbr)
                            visited[nbr[0], nbr[1]] = 1
                last_neighbour = pixel
            if not not_edge :
                edges.append(edge)
        ("computed the edges for junction {}".format(i-1))
        edges_list.append(edges)
        list = []
        for edge in edges :
            list.append(junction_coordinates.index(edge[-1]))
        junction_neighbours.append(list)
    return junction_neighbours, edges_list

            
def build_graph(skeleton, junction_coordinates : list):
    """
    Builds a graph with given junctions and following the edges determined by the skeleton (white pixels in `skeleton`).
    Parameters : 
    ---------------
    `skeleton` : image in black and white, [0, 255]
    `junction_coordinates` : list containing tuples of junctions coordinates in `skeleton`.

    Returns : 
    ---------------
    `G`: graph where the value of a vertex corresponds to the index of corresponding junction in `junction_coordinates`.
    """
    G=nx.Graph()  
    width = len(skeleton[0])
    height = len(skeleton)
    junction_indices = [junction[0]*width + junction[1] for junction in junction_coordinates]
    for i in range(height) :
        for j in range(width) :
            if skeleton[i, j] == 255 : 
                G.add_node(i*width+j)
                n1 = [(i, j-1), (i-1, j), (i+1, j), (i, j+1)]
                n2 = [(i-1, j-1), (i+1, j-1), (i-1, j+1), (i+1,j+1)]
                for n in n1 : 
                    if skeleton[n[0], n[1]] == 255 : 
                        G.add_edge(i*width+j, n[0]*width+n[1], weight = 1)
                for n in n2 : 
                    if skeleton[n[0], n[1]] == 255 : 
                        #diagonale must be penalizsd in order to favor neighbours in 4-connexity
                        G.add_edge(i*width+j, n[0]*width+n[1], weight = 3)
    return G, junction_indices

def get_edges(G, junction_indices) :
    edges_list = []
    for i in range(len(junction_indices)-1):
        print("junction {}".format(i))
        edges = []
        source = junction_indices[i]
        for j in range(i+1, len(junction_indices)):
            target = junction_indices[j]
            if nx.has_path(G, source, target) :
                #paths = nx.all_simple_paths(G, source, target, 400)
                print("found all paths :  ")
                i = 0
                print("here and")
                for path in nx.all_simple_paths(G, source, target, 200) : 
                    print("there")
                    print(i)
                    i += 1
                    if set(set(junction_indices) & set(path)) == set([source, target]) : 
                        #true edge : it doesn't pass by another junction
                        print("adding...")
                        edges.append(path)
                        print("added!")
        edges_list.append(edges)
    return edges_list



def local_extrema(edge : list, start : int, end : int, width : int) :
    max_dist = 0
    extrema = (start//width, start%width)
    ux = end//width - start//width
    uy = end%width - start%width
    norm_u = np.sqrt(ux**2 + uy**2)
    if norm_u > 0 :
        for point in edge : 
            vx = point//width - start//width
            vy = point%width - start%width
            dist = (ux*vy - uy*vx)/norm_u
            if dist > max_dist : 
                extrema = (point//width, point%width)
                max_dist = dist
    return extrema


def get_control_points_bezier(edge : list, P0 : int, P3 : int, width : int) -> list:
    """
    Computes the control points for a Bezier curve by using start and end points and furthest point from segment[start, end].
    
    Parameters : 
    ----------------
    `edge` : list of coordinates of pixels belonging to the edge
    `P0` : start point
    `P3` : end point
    
    Returns : 
    ----------------
    `coordinates` : list of coordinates of P0, P1, P2, P3 in skeleton. P0 and P3 have integer coordinates, P1 and P2 may not.
    """
    coordinates = [(P0//width, P0%width)]

    smoothness = 0.3
    x0 = P0//width
    y0 = P0%width
    x3 = P3//width
    y3 = P3%width
    extrema = local_extrema(edge, P0, P3, width)
    print(P0, P3, extrema)
    xm = extrema[0]
    ym = extrema[1]

    d01 = np.sqrt((xm-x0)**2 + (ym-y0)**2)
    d12 = np.sqrt((x3-xm)**2 + (y3-ym)**2)
    fa = smoothness*d01/(d01+d12)
    fb = smoothness*d12/(d01+d12) 
    p1x = xm - fa*(x3-x0) 
    p1y = ym - fa*(y3-y0)  
    p2x = xm + fb*(x3-x0)
    p2y = ym + fb*(y3-y0)  
    coordinates.append((p1x, p1y))
    coordinates.append((p2x, p2y))
    coordinates.append((x3, y3))
    return coordinates


def bezier_curve(tp : tuple, control_points : list):
    """
    Computes the Bezier curve determined by `control_points` at normalized position tp.
    
    Parameters : 
    ------------------
    `tp` : tuple of 2 elements containing normalized coordinates along initial curve.
    `control_points` : list of control points for the curve.
    
    Returns : 
    ------------------
    Coordinates of the point on the Bezier curve for normalized input tp."""

    p0 = control_points[0]
    p1 = control_points[1]
    p2 = control_points[2]
    p3 = control_points[3]
    x = (1-tp)**3*p0[0] + 3*tp*(1-tp)**2*p1[0] + 3*tp**2*(1-tp)*p2[0] + tp**3*p3[0]
    y = (1-tp)**3*p0[1] + 3*tp*(1-tp)**2*p1[1] + 3*tp**2*(1-tp)*p2[1] + tp**3*p3[1]
    return (x,y)


def fitting_error(edge : list, control_points : list, width) -> float:
    """
    Computes the fitting error of the Bezier curve for the edge.
    
    Parameters : 
    ------------------
    `edge`: list of coordinates of the pixels belonging to pixel. 
    `edge[0]` should be `control_points[0]` and `edge[-1]` should be equal to `control_points[3]`
    `control_points` : coordinates of 4 control points of the Bezier curve.
    
    Returns : 
    ------------------
    `fitting_error` : float 
    """
    w = 1 # for the moment I consider that thickness is invariant :
    #TODO : take into account line thickness
    sum = 0
    norm = len(edge)
    for i in range(norm) : 
        p = edge[i]
        p = (p//width, p%width)
        tp = i/norm
        B = bezier_curve(tp, control_points)
        sum += 0.5*((B[0] - p[0])**2 + (B[1] - p[1])**2)
        #print("norm {}, p {}, tp {}, B {}, sum {}".format(norm, p, tp, B, sum))
    return sum


def fit_bezier_curve(edge : list, start :int, end : int, junction_coordinates :list, width) :
    """
    
    Parameters : 
    ---------------------
    `edge` : all the coordinates of the pixels on the edge between `start` and `end`, in topological order.
    `start` : coordinates in image of starting point of `edge`
    `end` : coordinates in image of end point of `edge`
    `junction_coordinates` : list of the coordinates of all the junctions of the skeleton
    
    
    Returns : 
    --------------------
    `junction_coordinates` : the updated junction coordinates such that the fitting error of all Bezier curve is less than 2 pixels."""
    if len(edge)> 3 :     
        control_points = get_control_points_bezier(edge, start, end, width)
        err = fitting_error(edge, control_points, width)
        #print(err)
        if err > 20 :
            min_err = err
            split_index = 0
            for i in range(2,len(edge)-2) :
                err1 = fitting_error(edge[:i+1], get_control_points_bezier(edge[:i+1], edge[0], edge[i], width), width)
                err2 = fitting_error(edge[i:], get_control_points_bezier(edge[i:], edge[i], edge[-1], width), width)
                print("err1 {}, err {}".format(err1 + err2, min_err))
                if err1 + err2 < min_err :
                    split_index = i
                    min_err = err1 + err2
            split_point = edge[split_index]
            print(split_point)
            junction_coordinates.append((split_point//width, split_point%width))
            fit_bezier_curve(edge[:split_index+1], start, split_point, junction_coordinates, width)
            fit_bezier_curve(edge[split_index:], split_point, end, junction_coordinates, width)
    return junction_coordinates


def topological_graph(skeleton) -> tuple[list, list] : 
    """
    
    Parameters : 
    -----------------
    `skeleton` : image in black and white, [0, 255]
    
    Returns : 
    -----------------
    `nodes_coordinates` : the coordinates in `skeleton` of the nodes
    `nodes_neighbours` : same shape as `nodes_coordinates`. 
    ith element contains the list of nodes connected by an edge to pixel at nodes_coordinates[i].
    Neighbours correspond to their index in `nodes_coordinates`.
    `edges_list` : a list of same shape as `nodes_coordinates`, where edges_list[i] is the list of edges(all the pixels) for pixel at nodes_coordinates[i]. 
    Edges contain the coordinates of the pixels.

    """

    width = len(skeleton[0])
    print("getting junctions...")
    nodes_coordinates = get_junctions(skeleton)
    new_img = np.zeros(skeleton.shape)
    for node in nodes_coordinates :
        new_img[node[0], node[1]] = 255
    plt.imshow(new_img)
    plt.show()
    print("Done. Number of junctions identified : {} \n".format(len(nodes_coordinates)))
    print("getting neighbours and edges...")
    #junction_neighbours, edges_list = get_neighbours_and_edges(skeleton, nodes_coordinates)
    G, junctions_indices = build_graph(skeleton, nodes_coordinates)
    edges_list = get_edges(G, junctions_indices)
    print("Done.")
    new_nodes_coordinates = nodes_coordinates.copy()
    n = len(nodes_coordinates)
    for i in range(len(edges_list)) :
        start = junctions_indices[i]
        print("fitting edge {} out of {}".format(i, n))
        edges = edges_list[i]
        for edge in edges : 
            print("edge size : ", len(edge))
            new_nodes_coordinates = fit_bezier_curve(edge, start, edge[-1], new_nodes_coordinates, width)
    #nodes_neighbours, edges_list = get_neighbours_and_edges(skeleton, new_nodes_coordinates)
    G, junctions_indices = build_graph(skeleton, new_nodes_coordinates)
    print("Final number : {}".format(len(new_nodes_coordinates)))
    junctions_coordinates = [(i//width, i%width) for i in junctions_indices]
    return junctions_coordinates


def test() :
    name = 'boule'
    image = cv.imread("results/" + name +"_no_contour.png")
    image = cv.cvtColor(image, cv.COLOR_BGR2GRAY)

    skeleton = get_skeleton(image)
    plt.imshow(skeleton, cmap='gray')
    plt.show()

    nodes_coordinates = topological_graph(skeleton)
    new_img = np.zeros(skeleton.shape)
    for node in nodes_coordinates :
        new_img[node[0], node[1]] = 255
    plt.imshow(new_img)
    plt.show()

test()






    
            
