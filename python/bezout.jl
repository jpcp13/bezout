function bezout_gcd(a, b)
	s, t, r = 0, 1, b
	old_s, old_t, old_r = 1, 0, a
	while (r != 0)
	    quotient = div(old_r, r)
	    old_s, s = s, old_s - quotient*s
		old_t, t = t, old_t - quotient*t
		old_r, r = r, old_r - quotient*r
	end
    return 	old_s, old_t, old_r
end


function process_lines(i0, i1, j0)
	a = B[i0, j0]
	b = B[i1, j0]
	x, y, d = bezout_gcd(a, b)
	givens = reshape([x, -div(b, d), y, div(a, d)], (2, 2))
	if (d != 0)
		for k = 0:n
			B[[D*k + i0, D*k + i1], :] = givens*B[[D*k + i0, D*k + i1], :]
		end
	end
end

function process_block(block_start, block_end)
	for j0 = block_start:block_end
		for i1 = block_end:-1:j0+1
			print(i1)
			print(j0)
			process_lines(i1-1, i1, j0)
		end
	end
end


B = readdlm("bb.txt", BigInt)
D = size(B, 2)
n = div(size(B, 1), D) - 1
d = round(Int, (D/factorial(n))^(1/n))
block_size = div(D, n*d)


for l = 1:n*d
	block_end = block_size*l
	block_start = block_end - block_size + 1
	process_block(block_start, block_end)
end

writedlm("bb1_jl.txt", B)









