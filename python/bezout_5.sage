import numpy as np
import scipy.linalg as la
import time
import bezout_5 as bz

#TEX_DIR = '/home/jp/Documents/Bezout/bezout/tex/txt'
TEX_DIR = '../tex/txt'

deg = [1, 1, 1, 1, 1, 1]

with open(TEX_DIR+'/deg.txt', 'w') as f:
    f.write(str(deg))
    
t = 10
m = 16000
n = len(deg)

Field = QQ
R = PolynomialRing(Field, 'x', n)
x = R.gens()
xx = [x[0]**0] + list(x)

fshape = [d+1 for d in deg]
dx = [(i+1)*deg[i] for i in range(n)]
dy = [(n-i)*deg[i] for i in range(n)]
fn, Dx, Dy = factorial(n), prod(dx), prod(dy)
with open(TEX_DIR+'/Dx.txt', 'w') as f:
    f.write("{0:d}".format(Dx))


P = [bz.rand_poly(n-1, m, deg, t, x) for i in range(n)] + xx
#P = load('P.sobj')

save(P, 'P')
bz.P2txt(n, deg, P, TEX_DIR)
F = [bz.poly2prism(fshape, p) for p in P]


t = time.clock()
Gx, Gy, Hx, Hy = bz._GH(n, fn, deg, dx, dy)
H, K = bz._HK(n, Hx, Hy)
J = bz._J(Dx, Dy, F, n, fshape, dx, dy, Gx, Gy)
C = bz._C(n, Dx, Dy, J)
B = bz._B(n, Dx, Dy, H, K, C)
B = bz.block_triang(n, Dx, Dy, dx, dy, deg, B)
construction_B_time = time.clock() - t
with open(TEX_DIR+'/construction_B_time.txt', 'w') as f:
    f.write("{0:.4f}".format(construction_B_time))

BB = []
for k in range(n+1):
        Bk = matrix(Field, B[k])
        BB.append(Bk[:, :])
bz.bz2txt(n, TEX_DIR, BB)

"""
computing rank of B0
"""
t = time.clock()
r0 = rank(BB[0])
rank_B0_time = time.clock() - t

with open(TEX_DIR+'/rank_B0.txt', 'w') as f:
    f.write("{0:d}".format(r0))
with open(TEX_DIR+'/rank_B0_time.txt', 'w') as f:
    f.write("{0:.4f}".format(rank_B0_time))
#print("sage_rank = {0:d}".format(r0))

bb = []
for k in range(n+1):
    bb.append(np.array(BB[k], dtype=float))

b0 = bb[0]
numpy_rank = np.linalg.matrix_rank(b0)
#print("numpy_rank = {0:d}".format(numpy_rank))

"""
reduction of Bezoutian matrices
"""
epsi = abs(b0).max()/1e5

t = time.clock()
nb_relations = 1
while nb_relations > 0:
    bb, r0, nb_relations = bz.iteration(bb, r0, epsi)
    #print bb[0].shape, r0, nb_relations
#print "bbt"
bbt = []
for k in range(n+1):
    bbt.append(bb[k].T)
#bbt = normaliz(bbt, epsi)
nb_relations = 1
while nb_relations > 0:
    bbt, r0, nb_relations = bz.iteration(bbt, r0, epsi)
    #print bbt[0].shape, r0, nb_relations
reductions_time = time.clock() - t
with open(TEX_DIR+'/bezout_dim.txt', 'w') as f:
    f.write("{0:d}".format(r0))
with open(TEX_DIR+'/reductions_time.txt', 'w') as f:
    f.write("{0:.4f}".format(reductions_time))

"""
Numerical compation of the roots
"""
t = time.clock()
bbt_roots = bz.bbt2roots(bbt, epsi, r0)
compute_roots_time = time.clock() - t
with open(TEX_DIR+'/compute_roots_time.txt', 'w') as f:
    f.write("{0:.4f}".format(compute_roots_time))

test_bbtr = bz.roots_test(P, x, bbt_roots)
test_roots = np.sort(test_bbtr[:])

np.savetxt(TEX_DIR+'/test_roots.txt', test_roots,  fmt='%1.3e')

hist, bin_edges = np.histogram(np.log10(test_roots))

with open(TEX_DIR+'/histogram.txt', 'w') as f:
    for k in range(len(hist)):
        nb_roots = hist[k]
        left_bin, right_bin = bin_edges[k], bin_edges[k+1]
        if nb_roots > 0:
            if k < len(hist)-1:
                f.write("$[{0:2.1f}, {1:2.1f}]$ & ${2:d}$\\\\\n".format(left_bin, right_bin, nb_roots))
            else:
                f.write("$[{0:2.1f}, {1:2.1f}]$ & ${2:d}$\n".format(left_bin, right_bin, nb_roots))

t = time.clock()
GB, grobner_dim = bz.compute_grobner(R, P, n)
grobner_time = time.clock() - t
bz.gb2txt(TEX_DIR, GB)

with open(TEX_DIR+'/grobner_dim.txt', 'w') as f:
    f.write("{0}".format(grobner_dim))

with open(TEX_DIR+'/grobner_time.txt', 'w') as f:
    f.write("{0:.4f}".format(grobner_time))

with open(TEX_DIR+'/bezout_exact.txt', 'w') as f:
    X = matrix(ZZ, Dx, 0)
    X_ortho = X.kernel().basis_matrix().LLL().transpose()
    Y, Y_time = bz.BB2Y(BB, X_ortho)
    f.write("nb de Y-relations = {0:d}\n".format(Y.nrows()))
    
    K = bz.Y2K(BB[0], X_ortho, Y)
    X_dim = Dx - K.nrows()
    f.write("X_dim = {0:d}\n".format(X_dim))

    Y_ortho = Y.right_kernel_matrix().LLL()
    X, X_time = bz.BB2X(BB, Y_ortho)
    N = bz.X2N(BB[0], Y_ortho, X)
    X_ortho = X.kernel().basis_matrix().LLL().transpose()
    f.write("nb de X-relations = {0:d}\n".format(X.ncols()))
    
    XBBY = []
    for k in range(n+1):
        XBBY.append(X_ortho.transpose()*BB[k]*Y_ortho.transpose())
    bezout_exact_dim = rank(XBBY[0])
    f.write("bezout_exact_dim = {0:d}\n".format(bezout_exact_dim))

    K_ortho = XBBY[0].kernel().basis_matrix().right_kernel_matrix()
    N_ortho = N.kernel().basis_matrix().transpose()

    KBBN = []
    for k in range(n+1):
        KBBN.append(K_ortho*XBBY[k]*N_ortho)

    XX = []
    for k in range(n):
        xx = KBBN[0].solve_right(KBBN[k+1])
        XX.append(xx)

    test_XX = [norm(bz.X2p(XX, Field, p)) for p in P[:n]]
    
    f.write("test_XX = {0:s}".format(test_XX))

