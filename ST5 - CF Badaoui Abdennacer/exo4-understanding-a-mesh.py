# -*- coding: utf-8 -*-
"""
.. warning:: The explanations of the functions in this file and the details of
the programming methodology have been given during the lectures.
"""


# Python packages
import matplotlib.pyplot
import matplotlib.pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy
import os
import scipy.io
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import sys


# MRG packages
import solutions


def add_node_to_mesh(node_coords, p_elem2nodes, elem2nodes, nodeid_coords):
    # .. todo:: Modify the lines below to add one node to the mesh
    pass
    # return node_coords, p_elem2nodes, elem2nodes,


def add_elem_to_mesh(node_coords, p_elem2nodes, elem2nodes, elemid2nodes):
    # .. todo:: Modify the lines below to add one element to the mesh
    pass
    # return node_coords, p_elem2nodes, elem2nodes,


def remove_node_to_mesh(node_coords, p_elem2nodes, elem2nodes, nodeid):
    # .. todo:: Modify the lines below to remove one node to the mesh
    pass
    # return node_coords, p_elem2nodes, elem2nodes,


def remove_elem_to_mesh(node_coords, p_elem2nodes, elem2nodes, elemid):
    # .. todo:: Modify the lines below to remove one element to the mesh
    pass
    # return node_coords, p_elem2nodes, elem2nodes,


def compute_barycenter_of_element(node_coords, p_elem2nodes, elem2nodes):
    # ..todo: Modify the lines below to compute the barycenter of one element
    pass
    # return elem_coords


def compute_aspect_ratio_of_element(node_coords, p_elem2nodes, elem2nodes):
    # .. todo:: Modify the lines below to compute the quality criteria of all the elements
    pass
    # return elem_quality


def compute_edge_length_factor_of_element(node_coords, p_elem2nodes, elem2nodes):
    # .. todo:: Modify the lines below to compute the quality criteria of all the elements
    pass
    # return elem_quality


def compute_pointedness_of_element(node_coords, p_elem2nodes, elem2nodes):
    # .. todo:: Modify the lines below to compute the quality criteria of all the elements
    pass
    # return elem_quality


def run_exercise_a():
    """Generate grid with quadrangles.
    """
    # -- generate grid with quadrangles
    xmin, xmax, ymin, ymax = 0.0, 1.0, 0.0, 1.0
    nelemsx, nelemsy = 10, 10
    nelems = nelemsx * nelemsy
    nnodes = (nelemsx + 1) * (nelemsy + 1)
    # .. todo:: Modify the line below to call to generate a grid with quadrangles
    # p_elem2nodes, elem2nodes, node_coords, node_l2g = my_set_quadmesh(...)
    # .. note:: If you do not succeed, uncomment the following line to access the solution
    node_coords, node_l2g, p_elem2nodes, elem2nodes = solutions._set_quadmesh(xmin, xmax, ymin, ymax, nelemsx, nelemsy)

    # -- plot mesh
    fig = matplotlib.pyplot.figure(1)
    ax = matplotlib.pyplot.subplot(1, 1, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    solutions._plot_mesh(p_elem2nodes, elem2nodes, node_coords, color='yellow')
    node = 80
    solutions._plot_node(p_elem2nodes, elem2nodes, node_coords, node, color='red', marker='o')
    elem = 70
    solutions._plot_elem(p_elem2nodes, elem2nodes, node_coords, elem, color='orange')
    matplotlib.pyplot.show()

    return


def run_exercise_b():
    """Generate grid with triangles.
    """
    # -- generate grid with triangles
    xmin, xmax, ymin, ymax = 0.0, 1.0, 0.0, 1.0
    nelemsx, nelemsy = 10, 10
    nelems = nelemsx * nelemsy * 2
    nnodes = (nelemsx + 1) * (nelemsy + 1)
    node_coords, node_l2g, p_elem2nodes, elem2nodes = solutions._set_trimesh(xmin, xmax, ymin, ymax, nelemsx, nelemsy)
    nodes_on_boundary = solutions._set_trimesh_boundary(nelemsx, nelemsy)

    # -- plot mesh
    fig = matplotlib.pyplot.figure(1)
    ax = matplotlib.pyplot.subplot(1, 1, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    solutions._plot_mesh(p_elem2nodes, elem2nodes, node_coords, color='yellow')
    node = 80
    solutions._plot_node(p_elem2nodes, elem2nodes, node_coords, node, color='red', marker='o')
    elem = 70
    solutions._plot_elem(p_elem2nodes, elem2nodes, node_coords, elem, color='orange')
    matplotlib.pyplot.show()

    return


def run_exercise_c():
    pass


def run_exercise_d():
    pass


if __name__ == '__main__':

    run_exercise_a()
    run_exercise_b()
    run_exercise_c()
    run_exercise_d()
    print('End.')
