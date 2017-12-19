from sage.all import *
import scipy.fftpack as ft
import numpy as np
import scipy.linalg as la
import time
import os
import shutil
import matplotlib.pyplot as plt

def poly2prism(fshape, p):
	z = np.zeros(fshape)
	mns = p.monomials(); cfs = p.coefficients()
	for k in range(len(mns)):
		z[mns[k].degrees()] = cfs[k]
	return z

def _GH(n, fn, deg, dx, dy):
	Gx = []; Hx = []; Gy = []; Hy = []
	for j in range(n):
		dft = ft.fft(np.eye(2*fn*deg[j]))
		gx = dft[0::2*fn/(j+1), :deg[j]+1]
		gy = dft[1::2*fn/(n-j), :deg[j]+1]
		hx = dft[0::2*fn/(j+1), :dx[j]]
		hy = dft[1::2*fn/(n-j), :dy[j]]
		Gx.append(gx); Hx.append(hx); Gy.append(gy); Hy.append(hy)
	return Gx, Gy, Hx, Hy

def _HK(n, Hx, Hy):
	H = 1; K = 1
	for k in range(n):
		H = np.kron(Hx[k], H); K = np.kron(Hy[k], K)
	return H, K

def evaluate(F, n, fshape, dx, dy, Gx, Gy, i, j):
	f = F[i].transpose(range(n-1, -1, -1)).reshape(prod(fshape[j:]), prod(fshape[:j]))
	L = 1; R = 1
	for k in range(j):
		o = np.ones((dx[k], 1))
		L = np.kron(o, L); R = np.kron(Gy[k], R)
	for k in range(j, n):
		o = np.ones((dy[k], 1))
		R = np.kron(o, R); L = np.kron(Gx[k], L)
	return L.dot(f).dot(R.T)

def _J(Dx, Dy, F, n, fshape, dx, dy, Gx, Gy):
	J = np.empty((Dx, Dy, 2*n+1, n+1), dtype=complex)
	for i in range(2*n+1):
		for j in range(n+1):
			J[:, :, i, j] = evaluate(F, n, fshape, dx, dy, Gx, Gy, i, j)
	return J

def _C(n, Dx, Dy, J):
	C = np.empty((n+1, Dx, Dy), dtype=complex)
	indices_vol = range(n, 2*n+1)
	for i in range(n+1):
		indices = range(n) + [n+i]
		vol = np.linalg.det(J[:, :, indices_vol, :])
		num = np.linalg.det(J[:, :, indices, :])
		C[i] = num/vol
	return C

def _B(n, Dx, Dy, H, K, C):
	B = np.empty((n+1, Dx, Dy), dtype=complex)
	for i in range(n+1):
		HC  = np.conjugate(H.T).dot(C[i])/Dx
		KHC = np.conjugate(K.T).dot(HC.T)/Dy
		B[i] = KHC.T
	return np.around(B).real.astype(int)

def rand_poly(j, m, deg, t, x):
    # construit un polynome random en j+1 variables 0..j, de degres max(deg, m)
    p = 0
    for k in range(min(deg[j], m) + 1):
        if j > 0:
            p += rand_poly(j-1, m-k, deg, t, x)*x[j]**k
        else:
            coeff = ZZ.random_element(-t, t)*int(random() + 0.30)
            #coeff = Field.random_element()*int(random()*1.5)
            p += coeff*x[j]**k
    return p

def nzrows(m):
    nr = m.nrows()
    nzr = []
    for i in range(nr):
        if m[i, :].rank() > 0:
            nzr.append(i)
    return nzr

def mat2txt(mat, mat_name):
    with open(mat_name+".txt", "w") as f:
        mat_dict = mat.dict()
        for ij in mat_dict:
            c = mat_dict[ij]
            #~ if c != 1:
            f.write(str(c)+'\n')

def P2txt(n, deg, P, directory):
    deg_str = ''.join(str(e) for e in deg)
    local_dir = directory+'/P_'+deg_str
    if 'P_'+deg_str in os.listdir(directory):
        shutil.rmtree(local_dir)
    os.mkdir(local_dir)
    #~ pol2text(n, P, directory+'/P')
    #~ pol2tex(n, P, directory+'/P')
    for i in range(n):
        p = P[i]
        pdict = p.dict()
        with open(local_dir+'/p_'+str(i)+'.txt', 'w') as f:
            for key, val in pdict.items():
                f.write(str(key) + ', ' + str(val) + '\n')

def bz2txt(n, directory, BB):
    for k in range(n+1):
        m = BB[k]
        mat2txt(m, directory+'/BB/bb'+str(k))

def gb2txt(directory, GB):
	for k in range(len(GB)):
		p = GB[k]
		mat2txt(p, directory+'/GB/gb'+str(k))

def BB2roots(bb):
    qrBB = []
    n = len(bb)-1
    for k in range(n+1):
        qrBB.append(np.array(bb[k]))
    b0 = qrBB[0]
    bezout_dim = b0.shape[0]
    chow_mat = np.zeros((bezout_dim, bezout_dim))
    for k in range(n):
        chow_mat += np.random.randn()*qrBB[k+1]
    chowchow, b0b0, Q, Z = la.qz(chow_mat, b0, output='complex')
    qzBB = []
    for k in range(n+1):
        bk = qrBB[k]
        qzBB.append( Q.T.conjugate().dot(bk).dot(Z) )
    roots = np.zeros((bezout_dim, n), dtype=complex)
    for j in range(n):
        roots[:, j] = np.diag(qzBB[j+1])/np.diag(qzBB[0])
    return roots


def evalP(P, x, root):
	n = len(root)
	px = np.zeros((n), dtype='complex')
	Jx = np.zeros((n, n), dtype='complex')
	for i in range(n):
		pi = P[i]
		px[i] = pi(tuple(root))
		for j in range(n):
			djpi = P[i].derivative(x[j])
			Jx[i, j] = djpi(tuple(root))
	return px, Jx

def roots_test(P, x, roots):
	dim, n = roots.shape
	test = np.zeros((dim))
	test_r = np.zeros((dim))
	for k in range(dim):
		root = roots[k, :]
		px, Jx = evalP(P, x, root)
		iJpx = np.linalg.solve(Jx, px)
		test[k] = np.linalg.norm(iJpx)
		test_r[k] = np.linalg.norm(iJpx)/np.linalg.norm(root)
	return test_r

def compute_grobner(R, P, n):
    I = R.ideal(P[:n])
    GB = I.groebner_basis()
    grobner_dim = I.vector_space_dimension()
    return GB, grobner_dim

def X2m(XX, v, m):
    n = len(XX)
    r = v
    degs = m.degrees()
    for j in range(n):
        Xj = XX[j]
        for k in range(degs[j]):
            r *= Xj
    return r

def X2p(XX, Field, p):
    dim = XX[0].nrows()
    V = VectorSpace(Field , dim)
    pm = p.monomials()
    pc = p.coefficients()
    v = V.random_element()
    s = vector(Field, dim)
    for k in range(len(pc)):
        s += pc[k]*X2m(XX, v, pm[k])
    return s

def find_Y(BB, B0Y, Field, n):
    nr, nc = BB[0].nrows(), BB[0].ncols()
    while True:
        old_nbr = B0Y.nrows()
        print old_nbr
        K_lambda = B0Y.kernel().basis_matrix()
        K = K_lambda[:, 0:nc]
        Y = matrix(Field, 0, nc)
        for k in range(1, n+1):
            Yk = K*BB[k]
            Y = block_matrix(2, 1, [Y, Yk]).LLL()
            nzr = nzrows(Y) 
            Y = Y[nzr, :]
        B0Y = block_matrix(2, 1, [BB[0], Y])
        if B0Y.nrows() == old_nbr:
            break
            
    return Y

def P2field(p, Field):
    coeffs, monoms = p.coefficients(), p.monomials()
    return sum([Field(coeffs[i])*monoms[i] for i in range(len(coeffs))])
