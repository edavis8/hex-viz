#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:05:17 2021

@author: admin
"""


"""HEXAGON STUFF"""

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm
def make_hex_matrix(n=8):

    
    #Define a hexagon by rows
    rows = []
    m=n
    node = 0
    while m<=2*n-1:
        new_row = []
        for a in range(m):
            new_row.append(node)
            node+=1
        rows+=[new_row]
        m+=1
    m-=2
    while m>=n:
        new_row = []
        for a in range(m):
            new_row.append(node)
            node+=1
        rows+=[new_row]
        m-=1
     
    dim = node
    
    #fill our adjacency matrix using the hexagon rows
    matrix = [[0 for a in range(dim)] for b in range(dim)]    


    for i in range(len(rows)):
        for j in range(len(rows[i])):
            try:
                matrix[rows[i][j]][rows[i][j+1]] = 1
                matrix[rows[i][j+1]][rows[i][j]] = 1
            except:
                pass
            if i< n-1:
                matrix[rows[i][j]][rows[i+1][j]] = 1
                matrix[rows[i+1][j]][rows[i][j]] = 1
                
                matrix[rows[i][j]][rows[i+1][j+1]] = 1
                matrix[rows[i+1][j+1]][rows[i][j]] = 1
            else:
                try:
                    matrix[rows[i][j]][rows[i+1][j]] = 1
                    matrix[rows[i+1][j]][rows[i][j]] = 1
                except:
                    pass
                if j-1>=0:
                    try:
                        matrix[rows[i][j]][rows[i+1][j-1]] = 1
                        matrix[rows[i+1][j-1]][rows[i][j]] = 1
                    except:
                        pass
                    
    matrix = np.array(matrix)
    return matrix, rows


def plot_hex_matrix(matrix,rows, k_matrix):
    fig, ax = plt.subplots(figsize=(10,10))
    ax.axis('off')
    
    #start on the bottom row of the hexagon
    n = len(rows[0])
    plot_x = []
    plot_y = []
    y = 0
    
    #generate the cart coordinates for each node
    for row in rows:
        x = 1/2 - (len(row)/2)*(1/((2*n-1)))
        for node in row:
            plot_x+= [x]
            plot_y+=[y]
            x+=1/(2*n-1)
        y += 1/(2*n-1)

    #normalize the weights so each line is close to weight 1
    k_matrix_adj = (4/n)*(non_zero(k_matrix)/2)*k_matrix/np.sum(k_matrix)

    #plot each link individually
    for k in range(len(matrix)):
        for l in range(len(matrix)):
            if matrix[k,l] ==1:
                ax.plot([plot_x[k], plot_x[l]], [plot_y[k], plot_y[l]], color = 'gray', linewidth = k_matrix_adj[k,l])
    

###get our first batch fo conductances (random)
def init_conductances(matrix, gamma):
    #number of nodes
    dim = len(matrix)
    
    #get a random matrix and make sure it's symmetric
    random_matrix = np.random.rand(dim,dim)
    random_matrix = (random_matrix+random_matrix.T)/2
    
    k_matrix = random_matrix*matrix
    
    #normalize
    k_matrix = (k_matrix/(np.sum(k_matrix)))


    return k_matrix
    

 
def node_current_init(matrix, i_0):
    #all nodes are sinks except the source
    currents = np.zeros(len(matrix))
    i_k = -i_0/(len(matrix)-1)
    currents[0] = i_0
    for i in range(1,len(currents)):
        currents[i] = i_k
    return currents


#sum over rows for diagonal
def diagonal(k_matrix):
    diagonal_matrix = np.zeros((len(k_matrix), len(k_matrix)))
    for i in range(len(k_matrix)):
        diagonal_matrix[i,i] = sum(k_matrix[i])
    return diagonal_matrix


def solve_U(currents, k_matrix, diagonal_matrix):
    A = diagonal_matrix-k_matrix
    b= currents
    U = np.linalg.solve(A,b)
    
    return U


def link_currents(matrix, U, k_matrix, currents):
    
    I = np.zeros((len(matrix), len(matrix)))
    
    for k in range(len(matrix)):
        for l in range(len(matrix)):
            I[k,l] = k_matrix[k,l]*(U[l]-U[k])
    
    return I

#get new conductances for next t
def update_conductances(k_matrix, I, gamma):
    k_matrix = np.zeros((len(k_matrix), len(k_matrix)))
    big_gamma = 2*gamma/(gamma+1)
    for k in range(len(k_matrix)):
        for l in range(len(k_matrix)):
            k_matrix[k,l] = abs(I[k,l])**-(big_gamma-2)
    k_matrix = (k_matrix/(np.sum(k_matrix)))

    return k_matrix
    

    
#iterate and plot result
def main():
    matrix, rows  = make_hex_matrix(n)
    
    currents = node_current_init(matrix, i_0)
    
    k_matrix_init = init_conductances(matrix, gamma)
    
    k_matrix = k_matrix_init

    for i in tqdm(range(epochs)):
        try:
            diagonal_matrix = diagonal(k_matrix)
            
            U = solve_U(currents, k_matrix, diagonal_matrix)
            
            I = link_currents(matrix, U, k_matrix, currents)
            
            k_matrix = update_conductances(k_matrix, I, gamma)
        except:
            break
        
    plot_hex_matrix(matrix, rows, k_matrix)


#debugging stuff
def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)   
 
def non_zero(matrix):
    egh = [a for a in matrix.flatten() if a !=0]
    return len(egh)


#plot the lattice with appropriate line weights
def plot(gamma, n, i_0, epochs):
    matrix, rows  = make_hex_matrix(n)

    currents = node_current_init(matrix, i_0)

    k_matrix_init = init_conductances(matrix, gamma)

    k_matrix = k_matrix_init

    for i in tqdm(range(epochs)):
        try:
            diagonal_matrix = diagonal(k_matrix)

            U = solve_U(currents, k_matrix, diagonal_matrix)

            I = link_currents(matrix, U, k_matrix, currents)

            k_matrix = update_conductances(k_matrix, I, gamma)
        except:
            break

    plot_hex_matrix(matrix, rows, k_matrix)




if __name__ == '__main__':
    #our gamma value
    gamma = 0.5
    #side length of hexagon
    n=15
    #current source
    i_0 = 1000
    #time to run
    epochs= 50
    main()