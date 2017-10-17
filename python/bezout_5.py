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

def build_ixy(k, n, deg):
	ix = ()
	iy = ()
	for j in range(n):
		if j < k:
			ix = ix + (slice(None),)
			iy = (slice(None),) + iy
		elif j == k:
			ix = ix + (slice(-deg[k], None),)
			iy = (slice(-deg[n-1-k], None),) + iy
		else:
			ix = ix + (slice(None, -deg[j]),)
			iy = (slice(None, -deg[n-1-j]),) + iy
	return ix, iy

def permut(n, Dx, Dy, dx, dy, deg):
	aax = np.arange(Dx, dtype=int).reshape(dx[::-1]).transpose(range(n-1,-1,-1))
	aay = np.arange(Dy, dtype=int).reshape(dy[::-1]).transpose(range(n-1,-1,-1))
	ax = np.zeros(0, dtype=int)
	ay = np.zeros(0, dtype=int)
	for k in range(n):
		ix, iy = build_ixy(k, n, deg)
		tx = (np.arange(n)*(n-1) + k) % n
		ty = (n-1 - tx) % n
		ay = np.concatenate(( ay, np.transpose(aay[iy], ty).reshape(Dy/n) ))
		ax = np.concatenate(( ax, np.transpose(aax[ix], tx).reshape(Dx/n) ))
	return ax, ay

def block_triang(n, Dx, Dy, dx, dy, deg, B):
	ax, ay = permut(n, Dx, Dy, dx, dy, deg)
	for k in range(n+1):
		B[k] = np.fliplr(B[k][np.ix_(ax,ay)])
	return B

def rand_poly(j, m, deg, t, x):
    # construit un polynome random en j+1 variables 0..j, de degres max(deg, m)
    p = 0
    for k in range(min(deg[j], m) + 1):
        if j > 0:
            p += rand_poly(j-1, m-k, deg, t, x)*x[j]**k
        else:
            coeff = ZZ.random_element(-t, t)*int(random()+0.25)
            #coeff = Field.random_element()*int(random()*1.5)
            p += coeff*x[j]**k
    return p

def h_gcd(B):
    C = copy(B)
    for i in range(B.nrows()):
        g = gcd(B[i, :].list())
        if g != 0:
            C[i, :] = B[i, :]/g
    return C

def nzrows(m):
    nr = m.nrows()
    nzr = []
    for i in range(nr):
        if m[i, :].rank() > 0:
            nzr.append(i)
    return nzr

def Y2K(B, X_ortho, Y):
    nr = B.nrows()
    BY = block_matrix(2, 1, [B, Y])
    K_lambda = BY.kernel().basis_matrix()
    Y_tilde = K_lambda[:, 0:nr]
    nc = X_ortho.ncols()
    XoYt = block_matrix(2, 1, [X_ortho.transpose(), Y_tilde])
    K_lambda = XoYt.kernel().basis_matrix()
    K = K_lambda[:, 0:nc]  
    return K

def X2N(B, Y_ortho, X):
    nc = B.ncols()
    BX = block_matrix(1, 2, [B, X])
    N_mu = BX.right_kernel_matrix().transpose()
    X_tilde = N_mu[0:nc, :]
    nr = Y_ortho.nrows()
    YoXt = block_matrix(1, 2, [Y_ortho.transpose(), X_tilde])
    N_mu = YoXt.right_kernel_matrix().transpose()
    N = N_mu[0:nr, :]
    return N

def K2Y(K, BB, X_ortho):
    n = len(BB) - 1
    nc = BB[0].ncols()
    Y = matrix(ZZ, 0, nc)
    kox = K*X_ortho.transpose()
    for k in range(1, n+1):
        Yk = kox*BB[k]
        Y = block_matrix(2, 1, [Y, Yk]).LLL()
        nzr = nzrows(Y)
        Y = Y[nzr, :]
    return Y

def N2X(N, BB, Y_ortho):
    n = len(BB) - 1
    nr = BB[0].nrows()
    X = matrix(ZZ, nr, 0)
    yon = Y_ortho.transpose()*N
    for k in range(1, n+1):
        Xk = BB[k]*yon
        X = block_matrix(1, 2, [X, Xk])
        Xt = X.transpose().LLL()
        nzr = nzrows(Xt)
        X = Xt[nzr, :].transpose()
    return X

def BB2Y(BB, X_ortho):
    n = len(BB) -1
    nc = BB[0].ncols()
    B = BB[0]
    Y = matrix(ZZ, 0, nc)
    r, rn = 0, 1
    t = time.clock()
    while r < rn:
        r = Y.nrows()
        K = Y2K(B, X_ortho, Y)
        new_Y = K2Y(K, BB, X_ortho).LLL()
        nzr = nzrows(new_Y)
        rn = len(nzr)
        Y = new_Y[nzr, :]
    Y_time = time.clock() - t
    return Y, Y_time

def BB2X(BB, Y_ortho):
    n = len(BB) -1
    nr = BB[0].nrows()
    B = BB[0]
    X = matrix(ZZ, nr, 0)
    r, rn = 0, 1
    t = time.clock()
    while r < rn:
        r = X.ncols()
        N = X2N(B, Y_ortho, X)
        new_X = N2X(N, BB, Y_ortho)
        new_Xt = new_X.transpose().LLL()
        nzr = nzrows(new_Xt)
        rn = len(nzr)
        X = new_Xt[nzr, :].transpose()
    X_time = time.clock() - t
    return X, X_time

def null_space(B, epsi, *args):
    nr, nc = B.shape
    q, r, p = la.qr(B, pivoting=True)
    if len(args) > 0:
        rb = args[0]
    else:
        d = np.diag(r)
        rb = sum(np.abs(d) > epsi)
        #print "relations", d
    null = q.T[rb:, :]
    null_ortho = q.T[:rb, :]
    return null, null_ortho

def dim_inter(b_rows, y_rows):
    y, d, z = la.svd(b_rows.dot(y_rows.T))
    return d

def relations(ker, bb, epsi):
    n = len(bb) -1
    nr, nc = bb[0].shape
    y = np.zeros((0, nc))
    for k in range(1, n+1):
        rel = ker.dot(bb[k])
        y = np.concatenate((y, rel), axis=0)
    rows_ortho, rows = null_space(y.T, epsi)
    return rows, rows_ortho

def iteration(bb, r0, epsi):
    n = len(bb) -1
    bb_out = []
    ker, ker_ortho = null_space(bb[0], epsi, r0)
    b0_ortho, b0_rows = null_space(bb[0].T, epsi, r0)
    y_rows, y_ortho = relations(ker, bb, epsi)
    nb_relations = y_rows.shape[0]
    if nb_relations > 0:
        y, d, z = la.svd(b0_ortho.dot(y_rows.T))
        #print "trigo", d
        deficit = sum(np.abs(d) < 1e-10)
        r0_out = r0 - deficit
        for k in range(n+1):
            bb_out.append(bb[k].dot(y_ortho.T))
    else:
        r0_out = r0
        for k in range(n+1):
            bb_out.append(bb[k])
    return bb_out, r0_out, nb_relations

def mat2txt(mat, mat_name):
    f = open(mat_name+".txt", "w")
    mat_dict = mat.dict()
    for ij in mat_dict:
        c = mat_dict[ij]
        #~ if c != 1:
        f.write(str(c)+'\n')
    f.close()

def pol2tex(n, P, P_filename):
    sp = str(P[:n])
    sp = sp.replace('[', '')
    sp = sp.replace(']', '')
    sp = sp.replace('x', "x_")
    sp = sp.replace('*', '')
    sp = sp.replace(',', ",\\\\")
    with open(P_filename+".txt", 'w') as f:
        f.write(sp)
	
def pol2text(n, P, P_filename):
    with open(P_filename+".text", 'w') as f:
        for i in range(n):
            sp = str(P[i])
            f.write(sp + '\n')

def P2txt(n, deg, P, directory):
    deg_str = ''.join(str(e) for e in deg)
    local_dir = directory+'/P_'+deg_str
    if 'P_'+deg_str in os.listdir(directory):
        shutil.rmtree(local_dir)
    os.mkdir(local_dir)	
    save(P, local_dir+'/P')
    pol2text(n, P, local_dir+'/P')
    pol2tex(n, P, local_dir+'/P')
    for i in range(n):
        p = P[i]
        pdict = p.dict()
        with open(local_dir+'/p_'+str(i)+'.txt', 'w') as f:
            for key, val in pdict.items():
                f.write(str(key) + ', ' + str(val) + '\n')

def bz2txt(n, directory, BB):
    #~ os.mkdir(directory+'/BB')
    for k in range(n+1):
        m = BB[k]
        mat2txt(m, directory+'/BB/bb'+str(k))
    #~ mat2txt(ker, directory+'/BB/ker')
    #~ mat2txt(ker_ortho, directory+'/BB/ker_ortho')
    #~ mat2txt(relations, directory+'/BB/relations')

def dims2txt(directory, BB, bezout_dim, groebner_dim, Dx, deg):
    with open(directory+'/rank_B0.txt', 'w') as f:
        f.write(str(rank(BB[0])))
    with open(directory+'/bezout_dim.txt', 'w') as f:
        f.write(str(bezout_dim))
    with open(directory+'/groebner_dim.txt', 'w') as f:
        f.write(str(groebner_dim))
    with open(directory+'/Dx.txt', 'w') as f:
        f.write(str(Dx))
    with open(directory+'/deg.txt', 'w') as f:
        f.write(str(deg))

def time2txt(directory, bezout_reductions_time, groebner_basis_time):
    with open(directory+'/bezout_time.txt', 'w') as f:
        f.write(str(bezout_reductions_time))
    with open(directory+'/bezout_time.txt', 'w') as f:
        f.write(str(groebner_basis_time)+'\n')

def gb2txt(directory, GB):
	#~ os.mkdir(directory+'/GB')
	for k in range(len(GB)):
		p = GB[k]
		mat2txt(p, directory+'/GB/gb'+str(k))

def BB2qr(n, relations, ker, ker_ortho, BB):
	np_rel = np.array(relations)
	nbr = np_rel.shape[0]
	if nbr > 0:
		qq, rr = np.linalg.qr(np_rel.T, mode='complete')
		wt = qq[:, nbr:]
	else:
		wt = np.eye(np_rel.shape[1])
	np_ker = np.array(ker)
	np_ko = np.array(ker_ortho)
	npBB = []
	for k in range(n+1):
		np_Bk = np.array(BB[k])
		ko_Bk = np_ko.dot(np_Bk)
		ker_Bk = np_ker.dot(np_Bk)
		npBB.append(ko_Bk.dot(wt))
	q, r = la.qr(npBB[0].T, mode='full')
	qrBB = []
	for k in range(n+1):
		qrBB.append(npBB[k].dot(q))
	return qrBB

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

def bbt2roots(bbt, epsi, r0):
    right_ker, right_ortho = null_space(bbt[0], epsi, r0)
    left_ker, left_ortho = null_space(bbt[0].T, epsi, r0)
    qrBB = []
    n = len(bbt)-1
    for k in range(n+1):
        qrbbk = right_ortho.dot(bbt[k]).dot(left_ortho.T)
        qrBB.append(qrbbk)
    b0 = qrBB[0]
    bezout_dim = b0.shape[0]
    chow_mat = np.zeros((r0, r0))
    for k in range(n):
        chow_mat += np.random.randn()*qrBB[k+1]
    chowchow, b0b0, Q, Z = la.qz(chow_mat, b0, output='complex')
    qzBB = []
    for k in range(n+1):
        bk = qrBB[k]
        qzBB.append( Q.T.conjugate().dot(bk).dot(Z) )
    roots = np.zeros((r0, n), dtype=complex)
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

def roots_plot(directory, test_r):
	hist, bins = np.histogram(np.log10(test_r), bins=30)
	width = 0.7 * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	plt.clf()
	plt.bar(center, hist, align='center', width=width)
	plt.savefig(directory+'/histo_roots.png')


def compute_groebner(R, P, n):
    I = R.ideal(P[:n])
    GB = I.groebner_basis()
    I = R.ideal(P[:n])
    t = time.clock()
    GB = I.groebner_basis()
    groebner_time = time.clock() - t
    groebner_dim = I.vector_space_dimension()
    return GB, groebner_time, groebner_dim

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

