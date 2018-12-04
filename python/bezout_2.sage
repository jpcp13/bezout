import numpy as np
import scipy.linalg as la
import time
import bezout_2 as bz

TEX_DIR = '../tex/txt'

deg = [4, 3, 5]
with open(TEX_DIR+'/deg.txt', 'w') as f:
    f.write(str(deg))
    
t = 5
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

P = load('P_'+''.join(str(e) for e in deg)+'.sobj')
#~ P = [bz.rand_poly(n-1, m, deg, t, x) for i in range(n)] + xx
#~ save(P, 'P_'+''.join(str(e) for e in deg))
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

BB = []
for k in range(n+1):
        Bk = matrix(Field, B[k])
        BB.append(Bk)
    
bezout_size = sum([float(os.path.getsize(TEX_DIR+'/BB/'+f)) \
                   for f in os.listdir(TEX_DIR+'/BB')])

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

t = time.clock()
BBN = bz.XY_reduct(BB, Field, n)
bezout_exact_dim = rank(BBN[0])
reductions_time = time.clock() - t

with open(TEX_DIR+'/reductions_time.txt', 'w') as f:
    f.write("{0:.4f}".format(reductions_time))
with open(TEX_DIR+'/bezout_exact_dim.txt', 'w') as f:
    f.write("{0:d}".format(bezout_exact_dim))

Field = GF(next_prime(2000))
test_XX = bz.Y_test(BBN, n, Field, P)
print("test_XX = {0:s}".format(test_XX))

"""
Numerical computation of the roots
"""
import numpy as np
bb = [np.array(BBN[k]) for k in range(n+1)]

t = time.clock()
roots = bz.qz_BB2roots(bb)
compute_roots_time = time.clock() - t

with open(TEX_DIR+'/compute_roots_time.txt', 'w') as f:
    f.write("{0:.4f}".format(compute_roots_time))

test_roots = bz.roots_test(P, x, roots)
hist, bin_edges = np.histogram(np.log10(test_roots + 1e-16), bins='scott')
with open(TEX_DIR+'/histogram.txt', 'w') as f:
    for k in range(len(hist)):
        nb_roots = hist[k]
        left_bin, right_bin = bin_edges[k], bin_edges[k+1 ]
        if nb_roots > 0 :
            if k < len(hist)-1 :
                f.write("$[{0:2.1f}, {1:2.1f}]$ & ${2:d}$\\\\\n".format(left_bin, right_bin, nb_roots))
            else:
                f.write("$[{0:2.1f}, {1:2.1f}]$ & ${2:d}$\n".format(left_bin, right_bin, nb_roots))
hist, bin_edges

"""
Grobner computations 
"""
t = time.clock()
GB, grobner_dim = bz.compute_grobner(R, P, n)
print("grobner_dim = {0}".format(grobner_dim))
grobner_time = time.clock() - t
bz.gb2txt(TEX_DIR, GB)
grobner_size = sum([float(os.path.getsize(TEX_DIR+'/GB/'+f)) for f in os.listdir(TEX_DIR+'/GB')])
with open(TEX_DIR+'/grobner_dim.txt', 'w') as f:
    f.write("{0}".format(grobner_dim))
with open(TEX_DIR+'/grobner_time.txt', 'w') as f:
    f.write("{0:.4f}".format(grobner_time))
with open(TEX_DIR+'/grobner_size.txt', 'w') as f:
    f.write("{0:.4f}".format(grobner_size/1000000))
    

