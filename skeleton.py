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
from itertools import islice

#### debugging variables
debugging_bezier = False
debugging_final_edges = False

def get_skeleton(img) :
    """
    This function takes as input the trapped ball segmented image. It then computes the skeleton by intersecting iteratively all the zones.
    Parameters : 
    -----------------
    `img` : the image for the skeleton extraction. Must be black lines on white background.
    
    Returns : 
    -----------------
    `skeleton` : same shape as `img`, with values in {0, 255}. Skeleton pixels are 255, background is 0.
    """
    classes = list(set(img.flatten().tolist())) #gets a list of the different zones identified
    print("There are {} classes".format(len(classes)))
    masks = []
    print("Computing masks...")
    for nr_class in classes :
        mask = 1*(img==nr_class)
        masks.append(mask)
    print("Masks computed.")
    superposition = np.zeros(img.shape)
    print("Computing intersections and skeleton...")
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
    Checks whether a certain pixel is a junction node. A pixel is considered junction if it has at least three neighbour in 4-connexity.
    
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
    Goes through `skeleton` and determines which pixels are junction pixels.
    Parameters : 
    --------------
    `skeleton` : image in black and white, [0, 255]

    Returns :
    --------------
    `junction_coordinates` : the coordinates of the junctions in the image (i,j).
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


def neighbours_4_connexity(img, p_i : int, p_j : int):
    """
    Computes the neighbours in 4-connexity in the image. 
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
    candidates = [(p_i-1, p_j), (p_i, p_j-1), (p_i, p_j+1), (p_i+1, p_j)]#, (p_i-1, p_j-1), (p_i-1, p_j+1), (p_i+1, p_j-1), (p_i+1, p_j+1)]
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
    width = len(skeleton[0])
    i = 1
    for junction in junction_coordinates : 
        print("junction {} out of {}".format(i, len(junction_coordinates)))
        i += 1
        #get neighbours : 
        edges = []
        junction_neighbours = neighbours_4_connexity(skeleton, junction[0], junction[1])
        nb_potential_edges = len(junction_neighbours)
        j = 1
        print("{} potential edges to check.".format(nb_potential_edges))
        for n in range(nb_potential_edges) : 
            edge = [junction[0]*width +junction[1], junction_neighbours[n][0]*width + junction_neighbours[n][1]]

            j += 1
            not_edge = False
            last_neighbour = (edge[-2]//width, edge[-2]%width)
            while (edge[-1]//width, edge[-1]%width) not in junction_coordinates and not_edge == False :
                #find neighbour of latest pixel of edge
                
                pixel = edge[-1]
                #print("edge length : ", len(edge))
                neighbours = neighbours_4_connexity(skeleton, pixel//width, pixel%width) # one of them should be edge[-2]
                if len(neighbours) < 2 : 
                    #must be an artifact
                    print("{} is in edge but has {} neighbours".format(pixel, len(neighbours)))
                    print(skeleton[pixel//width-1:pixel//width+2, pixel%width-1:pixel%width+2])
                    not_edge = True
                else : 
                    #print(neighbours, pixel, "last", last_neighbour)
                    #print(skeleton[pixel//width-1:pixel//width+2, pixel%width-1:pixel%width+2])
                    for nbr in neighbours :
                        if nbr != last_neighbour:
                            #print(nbr, edge[-2])
                            edge.append(nbr[0]*width + nbr[1])
                last_neighbour = (pixel//width, pixel%width)
            if not not_edge :
                edges.append(edge)
        #print("computed the edges for junction {}".format(i-1))
        edges_list.append(edges)
    return edges_list

            
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
    `junction_indices` : the indices of the junctions. Corresponds to the flattened coordinate (y*width + x).
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
                #n2 = [(i-1, j-1), (i+1, j-1), (i-1, j+1), (i+1,j+1)]
                for n in n1 : 
                    if skeleton[n[0], n[1]] == 255 : 
                        G.add_edge(i*width+j, n[0]*width+n[1], weight = 1000000)
                        #G.add_edge(n[0]*width+n[1],i*width+j,  weight = 10001)
                #for n in n2 : 
                 #   if skeleton[n[0], n[1]] == 255 : 
                  #      pass
                        #diagonale must be penalised in order to favor neighbours in 4-connexity
                        #G.add_edge(i*width+j, n[0]*width+n[1], weight = 10003)
                        #G.add_edge(n[0]*width+n[1],i*width+j,  weight = 10003)
    return G, junction_indices

    
def get_edges(G, junction_indices) :
    """
    
    Returns :
    -----------------
    `edges_list` : a list where `edges_list[i]` corresponds to the list of edges starting from `junction_indices[i]`.
    """
    edges_list = []
    for i in range(len(junction_indices)):
        print("junction {}".format(i))
        edges = []
        source = junction_indices[i]

        for j in range(i+1, len(junction_indices)):
            target = junction_indices[j]
            if target != source : 
                #computing the path between source and target
                for neighbour in nx.neighbors(G, source) :
                    #print('neighbour : ', neighbour)

                    G.edges[source,neighbour]["weight"] = 1
                    #G.edges[neighbour,source]["weight"] -= 10000

                    if nx.has_path(G, source, target) :
                        #print("there is a path")
                        #path = nx.dijkstra_path(G, neighbour, target, weight="weight")
                        #if set(set(junction_indices) & set(path)) == set([target, source]) : 
                            #true edge : it doesn't pass by another junction
                        paths = nx.all_simple_paths(G,source=source, target=target, cutoff = 300)
                        for path in paths : 
                            if set(set(junction_indices) & set(path)) == set([target, source]) : 
                                #print("adding {} and {}...".format(source, target))
                                edges.append(path)
                                #print("added!")
 
                    G.edges[source,neighbour]["weight"] = 1000000
                    #G.edges[neighbour, source]["weight"] += 10000
                #print("found {} edges".format(len(edges)))
        if len(edges) > 0 : 
            edges_list.append(edges)
    print("total : ", np.sum(len(edges_list[i]) for i in range(len(edges_list))))
    return edges_list



def local_extrema(edge : list, start : int, end : int, width : int) :
    """
    Takes the segment going from `start` and `end` and gets the point in image coordinates (i,j) """
    max_dist = 0
    extrema = (start//width, start%width)
    y_extr = start//width
    x_extr = start%width
    uy = abs(end//width - start//width)
    ux = abs(end%width - start%width)
    norm_u = np.sqrt(ux**2 + uy**2)
    if norm_u > 0 :
        for point in edge : 
            vy = abs(point//width - start//width)
            vx = abs(point%width - start%width)
            dist = abs(ux*vy - uy*vx)/norm_u #getting the distance between a point and a line through cross product
            if (dist > max_dist) or (dist ==  max_dist and y_extr + x_extr > point//width + point%width): 
                y_extr = point//width
                x_extr = point%width
                extrema = (y_extr, x_extr)
                max_dist = dist

    #print("local extrema : ", extrema)
    return extrema


def get_control_points_bezier(edge : list, P0 : int, P3 : int, width : int) -> list:
    """
    Computes the control points for a Bezier curve by using start (`P0`) and end (`P3`) points and furthest point from segment[`P0`, `P3`].
    
    Parameters : 
    ----------------
    `edge` : list of coordinates of pixels belonging to the edge
    `P0` : start point coordinates in flattened image
    `P3` : end point coordinates in flattened image
    
    Returns : 
    ----------------
    `coordinates` : list of coordinates of P0, P1, P2, P3 in skeleton. P0 and P3 have integer coordinates, P1 and P2 may not.
    """
    coordinates = [(P0//width, P0%width)]

    smoothness = 0.3
    y0 = P0//width
    x0 = P0%width
    y3 = P3//width
    x3 = P3%width
    if len(edge) > 5:
        extrema = local_extrema(edge, P0, P3, width)
        #print("P0, P3, extrema", P0, P3, extrema)
        ym = extrema[0]
        xm = extrema[1]
        if np.sqrt(x3**2+y3**2) < np.sqrt(x0**2 + y0**2) : 
        #tackles the issue of the edges being oriented, 
        # always compute the bezier curve by considering the vertex closest to origin as the starting point
            d12 = np.sqrt((ym-y0)**2 + (xm-x0)**2)
            d01 = np.sqrt((y3-ym)**2 + (x3-xm)**2)
        else : 
            d01 = np.sqrt((ym-y0)**2 + (xm-x0)**2)
            d12 = np.sqrt((y3-ym)**2 + (x3-xm)**2)
        fa = smoothness*d01/(d01+d12)
        fb = smoothness*d12/(d01+d12) 
        p1x = xm - fa*abs(x3-x0) 
        p1y = ym - fa*abs(y3-y0)  
        p2x = xm + fb*abs(x3-x0)
        p2y = ym + fb*abs(y3-y0)  
        #print("d01 {} d12 {} fa {} fb {} p1x {} p1y {} p2x {} p2y {}".format(d01, d12, fa, fb, p1x, p1y, p2x, p2y))

        coordinates.append((p1y, p1x))
        coordinates.append((p2y, p2x))
    else : 
        coordinates.append((y0, x0))
        coordinates.append((y3,x3))
    coordinates.append((y3, x3))
    return coordinates


def bezier_curve(tp : tuple, control_points : list):
    """
    Computes the Bezier curve determined by `control_points` at normalized position `tp`.
    
    Parameters : 
    ------------------
    `tp` : tuple of 2 elements containing normalized coordinates along initial curve.
    `control_points` : list of control points for the curve.
    
    Returns : 
    ------------------
    Coordinates of the point on the Bezier curve for normalized input tp."""
    degree = len(control_points)-1
    if degree == 3 : 
        p0 = control_points[0]
        p1 = control_points[1]
        p2 = control_points[2]
        p3 = control_points[3]
        #since edges are oriented in my construction (exist 2 times) 
        # I must get the same bezier curve for both edges and that means choosing one configuration for the polynom
        # I decide to always take the one obtained with the edge whose starting point is the closest to the origin
        if np.sqrt(p0[0]**2 + p0[1]**2) < np.sqrt(p3[0]**2+p3[1]**2):
            x = (1-tp)**3*p0[0] + 3*tp*(1-tp)**2*p1[0] + 3*tp**2*(1-tp)*p2[0] + tp**3*p3[0]
            y = (1-tp)**3*p0[1] + 3*tp*(1-tp)**2*p1[1] + 3*tp**2*(1-tp)*p2[1] + tp**3*p3[1]
        else : 
            x = (1-tp)**3*p3[0] + 3*tp*(1-tp)**2*p1[0] + 3*tp**2*(1-tp)*p2[0] + tp**3*p0[0]
            y = (1-tp)**3*p3[1] + 3*tp*(1-tp)**2*p1[1] + 3*tp**2*(1-tp)*p2[1] + tp**3*p0[1]
    elif degree == 2 : 
        p0 = control_points[0]
        p1 = control_points[1]
        p2 = control_points[2]
        if np.sqrt(p0[0]**2 + p0[1]**2) < np.sqrt(p2[0]**2+p2[1]**2):
            x = (1-tp)**2*p0[0] + 2*tp*(1-tp)*p1[0] + tp**2*p2[0]
            y = (1-tp)**2*p0[1] + 2*tp*(1-tp)**p1[1] + tp**2*p2[1]
        else : 
            x = (1-tp)**2*p2[0] + 2*tp*(1-tp)*p1[0]+ tp**2*p0[0]
            y = (1-tp)**2*p2[1] + 2*tp*(1-tp)**p1[1] + tp**2*p0[1]
    elif degree ==1 
        
    return (x,y)


def compute_bezier_curve(edges_list, junctions_indices, width):
    """
    Parameters : 
    ----------------

    Returns : 
    ----------------
    `curve` : a list of length 2 where `curve[0]` are the x coordinates and `curve[1]` are the y coordinates.
    """
    curve = [[],[]]
    for P0 in range(len(junctions_indices)) :
        edges = edges_list[P0]
        #print(len(edges))
        for edge in edges :
            x = []
            y = []
            if len(edge)>2:
                control_points = get_control_points_bezier(edge, edge[0], edge[-1], width)
                norm = len(edge)
                for i in range(norm) :
                    tp =i/norm
                    bezier = bezier_curve(tp, control_points)
                    x.append(bezier[1])
                    y.append(bezier[0])
                    curve[0].append(bezier[0])
                    curve[1].append(bezier[1])
            if debugging_bezier == True :
                plt.plot(x, y, 'bo', ms = 1)
                plt.plot([control_points[0][1], control_points[3][1]], [control_points[0][0], control_points[3][0] ], 'ro', ms = 5)
                plt.plot([control_points[1][1], control_points[2][1]], [control_points[1][0], control_points[2][0] ], 'go', ms = 5)
                ax = plt.gca()
                ax.set_xlim(0,width)
                ax.set_ylim(width,0)
                ax.set_aspect('equal', adjustable='box')        
                plt.show()

    return curve




def fitting_error(edge : list, control_points : list, width : int) -> float:
    """
    Computes the fitting error of the Bezier curve for `edge`.
    
    Parameters : 
    ------------------
    `edge`: list of coordinates of the pixels belonging to pixel. 
    `edge[0]` should be `control_points[0]` and `edge[-1]` should be equal to `control_points[3]`
    `control_points` : coordinates of 4 control points of the Bezier curve.
    `width` : width of the image, used to compute coordinates in the image.
    
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


def fit_bezier_curve(edge : list, start :int, end : int, junction_coordinates :list, width : int, threshold = 300) -> list :
    """
    Takes `edge` and recursively fits Bezier curve on different segments by splitting in new segments, in order to get a fitting error below `threshold`.
    Parameters : 
    ---------------------
    `edge` : all the coordinates of the pixels on the edge between `start` and `end`, in topological order.
    `start` : coordinates in image of starting point of `edge`
    `end` : coordinates in image of end point of `edge`
    `junction_coordinates` : list of the coordinates of all the junctions of the skeleton
    `width` : width of the image. Used to normalize.
    `threshold` : the error threshold. 
    
    Returns : 
    --------------------
    `junction_coordinates` : the updated junction coordinates such that the fitting error of all Bezier curve is less than 2 pixels.
    """

    if len(edge)> 10 :     
        control_points = get_control_points_bezier(edge, start, end, width)
        err = fitting_error(edge, control_points, width)
        #print(err)
        if err > threshold :
            min_err = err
            split_index = 0
            for i in range(2,len(edge)-2) :
                err1 = fitting_error(edge[:i+1], get_control_points_bezier(edge[:i+1], edge[0], edge[i], width), width)
                err2 = fitting_error(edge[i:], get_control_points_bezier(edge[i:], edge[i], edge[-1], width), width)
                #print("err1 {}, err {}".format(err1 + err2, min_err))
                if err1 + err2 < min_err :
                    split_index = i
                    min_err = err1 + err2
            split_point = edge[split_index]
            #print(split_point)
            junction_coordinates.append((split_point//width, split_point%width))
            if split_point != start and split_point != end : 
                fit_bezier_curve(edge[:split_index+1], start, split_point, junction_coordinates, width)
                fit_bezier_curve(edge[split_index:], split_point, end, junction_coordinates, width)
    return junction_coordinates


def topological_graph(skeleton, name) -> tuple[list, list] : 
    """
    
    Parameters : 
    -----------------
    `skeleton` : image in black and white, [0, 255]
    
    Returns : 
    -----------------
    `junction_coordinates` : the coordinates in `skeleton` of the nodes
    `nodes_neighbours` : same shape as `nodes_coordinates`. 
    ith element contains the list of nodes connected by an edge to pixel at nodes_coordinates[i].
    Neighbours correspond to their index in `nodes_coordinates`.
    `edges_list` : a list of same shape as `junction_coordinates`, where edges_list[i] is the list of edges(all the pixels) for pixel at `junction_coordinates[i]`. 
    Edges contain the coordinates of the pixels.

    """

    width = len(skeleton[0])
    print("getting junctions...")
    nodes_coordinates = get_junctions(skeleton)
    new_img = np.zeros(skeleton.shape)
    for node in nodes_coordinates :
        new_img[node[0], node[1]] = 255
    plt.imshow(new_img, cmap="gray")
    plt.title("junctions : {} identified".format(len(nodes_coordinates)))
    plt.show()
    cv.imwrite("results/"+name+"_initial_junctions.png", new_img)

    print("Done. Number of junctions identified : {} \n".format(len(nodes_coordinates)))

    print("getting neighbours and edges...")
    #junction_neighbours, edges_list = get_neighbours_and_edges(skeleton, nodes_coordinates)
    G, junctions_indices = build_graph(skeleton, nodes_coordinates)
    #edges_list = get_edges(G, junctions_indices)
    edges_list = get_neighbours_and_edges(skeleton, nodes_coordinates)
    #### Check number of final edges
    size = 0
    for edge_list in edges_list :
        size += len(edge_list)
    print("initial number of edges : ", size, " junctions : ", len(junctions_indices))
    new_img = np.zeros(skeleton.shape)
    for edges in edges_list :
        for edge in edges : 
            for pixel in edge : 
                new_img[pixel//width, pixel%width] = 255
    plt.imshow(new_img)
    plt.title("edges")
    plt.show()
    print("Done.")
    print("Going to fit Bezier curves on the edges to get the best nodes. \n")
    new_nodes_coordinates = nodes_coordinates.copy()
    n = len(nodes_coordinates)
    for i in range(len(edges_list)) :
        print("fitting edge {} out of {}".format(i, n))
        edges = edges_list[i]
        for edge in edges : 
            start = edge[0]
            #print("edge size : ", len(edge))
            new_nodes_coordinates = fit_bezier_curve(edge, start, edge[-1], new_nodes_coordinates, width)
    #nodes_neighbours, edges_list = get_neighbours_and_edges(skeleton, new_nodes_coordinates)
    G, junctions_indices = build_graph(skeleton, new_nodes_coordinates)
    print("Final number of nodes : {}".format(len(new_nodes_coordinates)))
    junctions_coordinates = [(i//width, i%width) for i in junctions_indices]
    final_edges = get_neighbours_and_edges(skeleton, junctions_coordinates)
    return junctions_coordinates, final_edges, junctions_indices


def test() :
    name = 'coffee_machine'
    image = cv.imread("results/" + name +"_no_contour.png")
    image = cv.cvtColor(image, cv.COLOR_BGR2GRAY)

    skeleton = get_skeleton(image)
    plt.imshow(skeleton, cmap='gray')
    plt.title("skeleton")
    plt.show()

    cv.imwrite("results/"+name+"_skeleton.png", skeleton)

    nodes_coordinates, edges_list, junctions_indices = topological_graph(skeleton, name)

#### Check number of final edges
    if debugging_final_edges : 
        nodes_coordinates = np.array(nodes_coordinates)
        size = 0
        max = 1
        nb = 1
        for edge_list in edges_list :
            size += len(edge_list)
            if size == max : 
                nb += 1
            if len(edge_list) > max : 
                nb = 1
                max = size
        print("final number of edges and junction and max nb neighbours: {} {} {} {}".format(size, len(junctions_indices), max, nb))
        fig = plt.figure()
        for edge_list in edges_list :
            for edge in edge_list :
                for point in edge : 
                    plt.plot(point%len(skeleton[0]), point//len(skeleton[0]), 'bo', ms = 1)
                
                plt.plot(nodes_coordinates[:,1], nodes_coordinates[:,0], 'ro', ms = 1)
                ax = plt.gca()
                ax.set_xlim(0,len(skeleton[0]))
                ax.set_ylim(len(skeleton),0)
                ax.set_aspect('equal', adjustable='box')
                ax.set_title("final edges")
                plt.show()


    new_img = np.zeros(skeleton.shape)
    for node in nodes_coordinates :
        new_img[node[0], node[1]] = 255
    plt.imshow(new_img)
    plt.title("final nodes")
    plt.show()
    cv.imwrite("results/"+name+"_nodes_bezier.png", new_img)
    curves = compute_bezier_curve(edges_list, junctions_indices, len(skeleton[0]))
    fig = plt.figure()
    plt.plot(curves[1], curves[0], 'bo', ms=0.08)
    ax = plt.gca()
    ax.set_xlim(0,len(skeleton[0]))
    ax.set_ylim(len(skeleton),0)
    ax.set_aspect('equal', adjustable='box')
    plt.title("bezier curves")
    plt.savefig("results/"+name+"_bezier_curves.png")
    plt.show()

test()


###debug control points and Bezier

def debug_bezier() :
    img = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for i in range(10)]
    img[1] = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
    img[2] = [1, 0, 1, 1, 0, 0, 0, 0, 0, 0]
    img[3] = [1, 0, 0, 1, 1, 0, 0, 0, 0, 0]
    img[4] = [1, 0, 0, 0, 1, 1, 0, 0, 0, 0]
    img[5] = [1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
    img[6] = [1, 0, 0, 0, 0, 0, 1, 1, 0, 0]
    plt.imshow(img)
    plt.show()

    edge = [61, 51, 41, 31, 21, 11, 12, 13, 23, 24, 34, 35, 45, 46, 56, 57, 67, 68]
    edges_list = [[[61, 51, 41, 31, 21, 11, 12, 13, 23, 24, 34, 35, 45, 46, 56, 57, 67, 68]]]
    edge.reverse()
    edges_list.append([edge])
    print(edges_list)
    points = get_control_points_bezier(edge, 68, 61, 10)
    print("points for reversed :", points)
    points2 = get_control_points_bezier(edges_list[0][0], 61, 68, 10)
    print("points : ", points2)
    curves = compute_bezier_curve(edges_list, [61, 68], 10)
    fig = plt.figure()
    n = len(curves[0])//2
    plt.plot(curves[1][0:n], curves[0][0:n], 'bo', ms=1)
    plt.plot(curves[1][n:], curves[0][n:], 'ro', ms=1)
    plt.plot(points[1], points[0], 'go', ms = 8)
    plt.plot(points2[1], points2[0], 'co', ms = 8)
    ax = plt.gca()
    ax.set_xlim(0,10)
    ax.set_ylim(10,0)
    ax.set_aspect('equal', adjustable='box')
    plt.title("bezier curves")
    plt.show()

#debug_bezier()



    
            
