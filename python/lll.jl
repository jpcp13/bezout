Pkg.add("LLLplus")
using LLLplus

# Time LLL, VBLAST decomposition of a complex matrix with randn entries 
N = 1000;
H = randn(N,N) + im*randn(N,N);
println("Testing LLL on $(N)x$(N) complex matrix...")
@time (B,T,Q,R) = lll(H);
M = 200;
println("Testing VBLAST on $(M)x$(M) chunk of same matrix...")
@time (W,P,B) = vblast(H[1:M,1:M]);

# Time LLL, Seysen decompositions of a 100x100 Int128 matrix with
# rand entries distributed uniformly between -100:100

BB = readdlm("bb.txt")
D = size(BB, 2)
H = BB[1:D, 1:D]


println("Testing VBLAST on $(M)x$(M) chunk of same matrix...")
@time (W,P,B) = vblast(H);
