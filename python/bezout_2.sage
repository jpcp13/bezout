import numpy as np
import scipy.linalg as la
import time
import bezout_2 as bz

#TEX_DIR = '/home/jp/Documents/Bezout/bezout/tex/txt'
TEX_DIR = '../tex/txt'

deg = [2, 2, 2, 2]

with open(TEX_DIR+'/deg.txt', 'w') as f:
    f.write(str(deg))
    
t = 12
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

#~ P = load('P_'+''.join(str(e) for e in deg)+'.sobj')
P = [bz.rand_poly(n-1, m, deg, t, x) for i in range(n)] + xx

save(P, 'P_'+''.join(str(e) for e in deg))
bz.P2txt(n, deg, P, TEX_DIR)
F = [bz.poly2prism(fshape, p) for p in P]

t = time.clock()
Gx, Gy, Hx, Hy = bz._GH(n, fn, deg, dx, dy)
H, K = bz._HK(n, Hx, Hy)
J = bz._J(Dx, Dy, F, n, fshape, dx, dy, Gx, Gy)
C = bz._C(n, Dx, Dy, J)
B = bz._B(n, Dx, Dy, H, K, C)

construction_B_time = time.clock() - t
with open(TEX_DIR+'/construction_B_time.txt', 'w') as f:
    f.write("{0:.4f}".format(construction_B_time))

#B = bz.block_triang(n, Dx, Dy, dx, dy, deg, B)
BB = []
for k in range(n+1):
        Bk = matrix(Field, B[k])
        BB.append(Bk[:, :])

bezout_size = sum([float(os.path.getsize(TEX_DIR+'/BB/'+f)) for f in os.listdir(TEX_DIR+'/BB')])
with open(TEX_DIR+'/bezout_size.txt', 'w') as f:
	f.write("{0:.4f}".format(float(bezout_size)/1000000))
	

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
print("sage_rank = {0:d}".format(r0))

bb = []
for k in range(n+1):
    bb.append(np.array(BB[k], dtype=float))

b0 = bb[0]
numpy_rank = np.linalg.matrix_rank(b0)
print("numpy_rank = {0:d}".format(numpy_rank))


"""
reduction process
"""
t = time.clock()
B0Y = BB[0]
Y = bz.find_Y(BB, B0Y, Field, n)
BBt = []
for k in range(n+1):
    BBt.append(BB[k].transpose())
B0Yt = B0Y.transpose()
X = bz.find_Y(BBt, B0Yt, Field, n)
print "-"*10
Y_ortho = Y.right_kernel_matrix()
bb0y = BB[0]*Y_ortho.transpose()
X_ortho = X.right_kernel_matrix()
xbb0 = X_ortho*BB[0]
xbb0y = X_ortho*BB[0]*Y_ortho.transpose()
bezout_exact_dim = rank(xbb0y)
print("bezout_exact_dim = {0:d}".format(bezout_exact_dim))
with open(TEX_DIR+'/bezout_dim.txt', 'w') as f:
    f.write("{0:d}".format(r0))
XBBY = []
for k in range(n+1):
    XBBY.append(X_ortho*BB[k]*Y_ortho.transpose())
K_ortho = XBBY[0].kernel().basis_matrix().right_kernel_matrix()
N_ortho = XBBY[0].right_kernel_matrix().right_kernel_matrix()
#N_ortho = N.kernel().basis_matrix().transpose()
KBBN = []
for k in range(n+1):
    KBBN.append(K_ortho*XBBY[k]*N_ortho.transpose())
bz.bz2txt(n, TEX_DIR, KBBN)
reductions_time = time.clock() - t
with open(TEX_DIR+'/reductions_time.txt', 'w') as f:
    f.write("{0:.4f}".format(reductions_time))


"""
Test using finite field arithmetic
"""
Field = GF(next_prime(200))
BBf = []
for k in range(n+1):
        Bk = matrix(Field, KBBN[k])
        BBf.append(Bk[:, :])
if rank(BBf[0]) == bezout_exact_dim:
    XX = []
    for k in range(n):
        xx = BBf[0].solve_right(BBf[k+1])
        XX.append(xx)
    Pf = [bz.P2field(p, Field) for p in P]
    test_XX = [bz.X2p(XX, Field, p) for p in Pf[:n]]
    #f.write("test_XX = {0:s}".format(test_XX))
    print("test_XX = {0:s}".format(test_XX))
else:
    print("rank deficiency !")


"""
Numerical compation of the roots
"""
t = time.clock()
bb = []
for k in range(n+1):
    bb.append(np.array(KBBN[k], dtype=float))
roots = bz.BB2roots(bb)
compute_roots_time = time.clock() - t
with open(TEX_DIR+'/compute_roots_time.txt', 'w') as f:
    f.write("{0:.4f}".format(compute_roots_time))
test_roots = bz.roots_test(P, x, roots)
print(test_roots)

hist, bin_edges = np.histogram(np.log10(test_roots), bins='scott')
with open(TEX_DIR+'/histogram.txt', 'w') as f:
    for k in range(len(hist)):
        nb_roots = hist[k]
        left_bin, right_bin = bin_edges[k], bin_edges[k+1 ]
        if nb_roots > 0 :
            if k < len(hist)-1 :
                f.write("$[{0:2.1f}, {1:2.1f}]$ & ${2:d}$\\\\\n".format(left_bin, right_bin, nb_roots))
            else:
                f.write("$[{0:2.1f}, {1:2.1f}]$ & ${2:d}$\n".format(left_bin, right_bin, nb_roots))


"""
Grobner computations 
"""
t = time.clock()
GB, grobner_dim = bz.compute_grobner(R, P, n)
print grobner_dim
grobner_time = time.clock() - t
bz.gb2txt(TEX_DIR, GB)
grobner_size = sum([float(os.path.getsize(TEX_DIR+'/GB/'+f)) for f in os.listdir(TEX_DIR+'/GB')])
with open(TEX_DIR+'/grobner_dim.txt', 'w') as f:
    f.write("{0}".format(grobner_dim))
with open(TEX_DIR+'/grobner_time.txt', 'w') as f:
    f.write("{0:.4f}".format(grobner_time))
with open(TEX_DIR+'/grobner_size.txt', 'w') as f:
    f.write("{0:.4f}".format(grobner_size/1000000))
