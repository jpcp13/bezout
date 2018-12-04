from sage.all import *
import scipy.fftpack as ft
import scipy.linalg as la
import numpy as np
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
    p = 0
    for k in range(min(deg[j], m) + 1):
        if j > 0:
            p += rand_poly(j-1, m-k, deg, t, x)*x[j]**k
        else:
            coeff = ZZ.random_element(-t, t)*int(random() + 0.2)
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

def pol2tex(n, P, P_filename):
    sp = str(P[:n])
    sp = sp.replace('[', '')
    sp = sp.replace(']', '')
    sp = sp.replace('x', "x_")
    sp = sp.replace('*', '')
    sp = sp.replace(',', ",\\\\")
    with open(P_filename+".txt", 'w') as f:
        f.write(sp)
        
def P2txt(n, deg, P, directory):
    deg_str = ''.join(str(e) for e in deg)
    local_dir = directory+'/P_'+deg_str
    if 'P_'+deg_str in os.listdir(directory):
        shutil.rmtree(local_dir)
    os.mkdir(local_dir)
    pol2tex(n, P, directory+'/P')
    for i in range(n):
        p = P[i]
        pdict = p.dict()
        with open(local_dir+'/p_'+str(i)+'.txt', 'w') as f:
            for key, val in pdict.items():
                f.write(str(key) + ', ' + str(val) + '\n')

def bz2txt(n, directory, BB):
    for f in os.listdir(directory+'/BB'):
		if f[-4:] == '.txt':
			os.remove(directory+'/BB/'+f)
    for k in range(n+1):
        m = BB[k]
        mat2txt(m, directory+'/BB/bb'+str(k))

def gb2txt(directory, GB):
	for f in os.listdir(directory+'/GB'):
		if f[-4:] == '.txt':
			os.remove(directory+'/GB/'+f)
	for k in range(len(GB)):
		p = GB[k]
		mat2txt(p, directory+'/GB/gb'+str(k))

def qz_BB2roots(bb):
    n = len(bb)-1
    bbt = [bb[k].T for k in range(n+1)] # transposition ! dim x DIM
    dim, DIM = bbt[0].shape
    q, r, p = la.qr(bbt[0], pivoting=True)
    qr_BB = []
    for k in range(n+1):
        qb = q.T.dot(bbt[k])
        qbp = qb[:, p]
        qr_BB.append(qbp[:, 0:dim])
    chow_mat = np.zeros((dim, dim))
    for k in range(n):
        chow_mat += np.random.randn()*qr_BB[k+1]
    b0 = qr_BB[0]
    chowchow, b0b0, Q, Z = la.qz(chow_mat, b0, output='complex')
    qzBB = []
    for k in range(n+1):
        qzBB.append( Q.T.conjugate().dot(qr_BB[k]).dot(Z) )
    roots = np.zeros((dim, n), dtype=complex)
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
	test = []
	test_r = np.zeros((dim))
	for k in range(dim):
		root = roots[k, :]
		px, Jx = evalP(P, x, root)
		try:
			iJpx = np.linalg.solve(Jx, px)
			test.append(np.linalg.norm(iJpx))
		except np.linalg.LinAlgError as err:
			if 'Singular matrix' in str(err):
				print("ooh")
			else:
				raise
	return np.array(test)
	
def Y_test(BBN, n, Field, P):
    def BBf(k):
        return matrix(Field, BBN[k])[:, :]
    if rank(BBf(0)) == BBN[0].ncols():
        XX = [BBf(0).solve_right(BBf(k+1)) for k in range(n)]
        Pf = [P2field(p, Field) for p in P]
        test_XX = [X2p(XX, Field, p).is_zero() for p in Pf[:n]]
        return test_XX
    else:
        print("rank deficiency ! Change finite field !")
        return None
        
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
        K = K_lambda[:, 0:nr]
        Y = matrix(Field, 0, nc)
        for k in range(1, n+1):
            Yk = K*BB[k]
            Y = block_matrix(2, 1, [Y, Yk]).row_space().basis_matrix()
            #~ nzr = nzrows(Y) 
            #~ Y = Y[nzr, :]
        B0Y = block_matrix(2, 1, [BB[0], Y])
        if B0Y.nrows() == old_nbr:
            break
    return Y

def Y_reduct(BB, Field, n):
    Y = find_Y(BB, BB[0], Field, n)
    E = Y.T.extended_echelon_form(subdivide=True)
    J = E[:, Y.nrows():]
    BJ = [BB[k]*J.T[:, Y.nrows():] for k in range(n+1)]
    K = BJ[0].kernel().basis_matrix()
    p = K.nonpivots()
    BJp = [BJ[k][p, :] for k in range(n+1)]
    return BJp

def XY_reduct(BB, Field, n):
    BJp = Y_reduct(BB, Field, n)
    BJpt = [b.T for b in BJp]
    BJpp = Y_reduct(BJpt, Field, n)
    BJppt = [b.T for b in BJpp]
    return BJppt

    
def P2field(p, Field):
    coeffs, monoms = p.coefficients(), p.monomials()
    return sum([Field(coeffs[i])*monoms[i] for i in range(len(coeffs))])
