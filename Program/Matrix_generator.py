import sympy as sp
import numpy as np
from scipy import sparse as sc


delta_tau = 0.02
#dim_y = 2


def create_a(lat_len, para1, para2):
    # Matrix Size
    J = para1
    Nt = int(para2)
    Ng = 2 * Nt * lat_len**2

    h = 1
    M = 10
    K = sp.atanh(sp.exp(-2 * delta_tau * h))

    # Create a map form x,y,tau,dir to o
    L = list()

    for dir in (0, 1):
        for t in range(int(Nt)):
            for y in range(lat_len):
                for x in range(lat_len):
                    L.append([x, y, t, dir])

    # Defining non-zero conditions
    def con_m(i,j):
        if L[i][0] == L[j][0] and L[i][1] == L[j][1] and L[i][2] == L[j][2] and L[i][3] == L[j][3]:
            return True

    def con_j(i,j):
        if L[i][3] == 0 and L[j][3] == 1 and L[i][2] == L[j][2] and ((
                                                                               L[i][0] == L[j][0] and (
                                                                               L[i][1] == L[j][1] or L[i][1] == (
                                                                               L[j][1] - 1) % lat_len)) or (
                                                                               L[i][0] == (L[j][0] + 1) % lat_len) and (
                                                                               L[i][1] == L[j][1] or L[i][1] == (
                                                                               L[j][1] - 1) % lat_len)):
            return True

        elif L[i][3] == 1 and L[j][3] == 0 and L[i][2] == L[j][2] and ((
                                                                               L[i][0] == L[j][0] and (
                                                                               L[i][1] == L[j][1] or L[i][1] == (
                                                                               L[j][1] + 1) % lat_len)) or (
                                                                               L[i][0] == (L[j][0] - 1) % lat_len) and (
                                                                               L[i][1] == L[j][1] or L[i][1] == (
                                                                               L[j][1] + 1) % lat_len)):
            return True

    def con_k(i,j):
        if L[i][0] == L[j][0] and L[i][1] == L[j][1] and L[i][3] == L[j][3] and (
                L[i][2] == (L[j][2] - 1) % Nt or L[i][2] == (L[j][2] + 1) % Nt):
            return True

    # Matrix gernerating function
    def generator(i, j):
        if con_m(i,j):
            return 2 * M
        elif con_j(i,j):
            return -delta_tau * J

        elif con_k(i,j):
            return K
        else:
            return 0

    A = sp.Matrix(Ng, Ng, generator)
    A = np.array(A).astype(np.float64)
    A_sparse = sc.csr_matrix(A)

    try:
        B = np.linalg.cholesky(A)
    except np.linalg.LinAlgError:
        print('not positiv definit')
    return A, A_sparse


def create_inv(lat_len, para1, para2):
    A, A_sparse = create_a(lat_len, para1, para2)
    matrix_inv = np.linalg.pinv(A)
    matrix_inv[abs(matrix_inv) < 1e-5] = 0

    return matrix_inv

