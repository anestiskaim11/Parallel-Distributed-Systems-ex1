### A Pluto.jl notebook ###
# v0.17.0

using Markdown
using InteractiveUtils

# ╔═╡ a15f08d6-8b69-4ed0-b22b-f920f98db371
begin
	using PlutoUI
	using SparseArrays
	using MatrixMarket
	using BenchmarkTools
end

# ╔═╡ 6339ca2a-fb15-41fe-b3e3-5e7fbbe5095f
#@author Anestis Kaimakamidis 9627 anestisk@ece.auth.gr
#Julia Implementation 
#NACA0015 com-Youtube dblp-2010 achieve better benchmarks than triangle counting on parallel implementations
#whereas mycielskian13 and belgium_osm are a bit slower than triangle counting function
#the code is the same as in the c files and explained in the report subitted

# ╔═╡ e58d2a23-799c-44fa-9956-9ecc27e4fc8f
begin
	A = MatrixMarket.mmread("NACA0015.mtx")
	A = sparse(A)
	view(A, :, :)
	#vals = findnz(A)
end

# ╔═╡ a9a46efc-aabf-4cb9-8980-cc74b245ace5
function triangle_counting(A)
	A .* (A * A)
end

# ╔═╡ 43b6eb91-e23c-4a2a-ab0a-5c219e9a965f
begin
	C = triangle_counting(A)
	e = ones( size(A,1) )
	c = C * e/2
	size(C)
end

# ╔═╡ 0b87b6c3-b074-4f7f-8016-725425091e7e
begin
	const rows = rowvals(A)
	const n, m = size(A)
	const vals = nonzeros(A)
	const range = []
	const thread_work = []
	const result_nzrange = zeros(n+1)
	const maxThreads = Threads.nthreads()
	for i in 1:n
		append!(range, nzrange(A, i).start)
	end
	append!(range, nzrange(A, n).stop+1)
	const nz = range[n+1]
	view(range, :)
end

# ╔═╡ dc12e43f-cb4c-465d-9da3-e4398ca041a4
function Do_Calc!(i, j, endrow, endcol)
	result = 0
	lasti = 0
	lastj = 0
	changed = false
	if i == endcol
		return result;
	end
	if j == endrow
		return result;
	end
	while(true)
		lasti = i
		lastj = j
		if rows[i] == rows[j]
			result += 1
			if j < endrow - 1
				j += 1
			end
			if i < endcol - 1
				i += 1
			end
		elseif rows[i] > rows[j]
			if j < endrow - 1
				j += 1
			elseif i < endcol - 1
				i += 1
			end
		elseif rows[i] < rows[j]
			if i < endcol - 1
				i += 1
			elseif j < endrow - 1
				j += 1
			end
		end
		if (i == endcol - 1 && j == endrow - 1 && (lasti < i || lastj < j))
			changed = true
		end
		if i >= endcol - 1 && j >= endrow - 1
			if rows[i] == rows[j] && changed
				result += 1
			end
			break
		end
	end	
	return result
end

# ╔═╡ c8e428df-a6ce-48de-a64a-ad0b0c9077d3
function sequential_masked_triangle_counting(A)
	result_rows = Int64[]
	result_nzrange = Int64[]
	result_nonzeros = Int64[]
	append!(result_nzrange, 1)
	result = 0
	counter = 0
	for i in 1:n
		for j in nzrange(A, i)
			result = Do_Calc!(range[i], range[rows[j]], range[rows[j]+1], range[i+1])
			if result != 0
				append!(result_nonzeros, result)
				append!(result_rows, rows[j])
				counter = counter + 1
			end
		end
		append!(result_nzrange , result_nzrange[i] + counter)
		counter = 0
	end
	
	return size(result_nonzeros)
end

# ╔═╡ b6c1cb40-6bb2-4fca-ba94-dd77986c8b13
function Threading!(i, res_nz, result_nonzeros, result_rows)
	for j in nzrange(A, i)
		result = Do_Calc!(range[i], range[rows[j]], range[rows[j]+1], range[i+1])
		if result != 0
    	 	res_nz[Threads.threadid()] += 1
		end
		result_nonzeros[j] = result
		result_rows[j] = rows[j]
	end
	result_nzrange[i] = range[i]
end

# ╔═╡ 8b6a63e7-5f92-4659-8dcb-239e29dc65b1
function threaded_masked_triangle_counting(A)
	result_rows = zeros(nz)
	result_nonzeros = zeros(nz)
	res_nz = zeros(maxThreads)
	results = zeros(maxThreads)

	Threads.@threads for i in 1:n
		Threading!(i, res_nz, result_nonzeros, result_rows)
	end
	return sum(res_nz)
end
		
	

# ╔═╡ fa133740-2aca-4948-87cc-45f2941e608f
begin
    @benchmark sequential_masked_triangle_counting($A)
end

# ╔═╡ 89a843ef-4eb6-4472-ab7a-f9b85e3c7200
begin
	@benchmark triangle_counting($A)
end

# ╔═╡ f2574652-9f34-4fb0-aed0-b5e09807f087
begin
	@benchmark threaded_masked_triangle_counting($A)
end

# ╔═╡ 5d781032-288d-4bfc-826b-31afdb8040ea
with_terminal() do
	sequential_masked_triangle_counting(A)
end

# ╔═╡ 702db5d8-7d85-41e6-ae14-1b3fbb84f7c9
with_terminal() do
	threaded_masked_triangle_counting(A)
end

# ╔═╡ a8bc2ff1-7ab0-4c6b-bf5b-31e43458ded0


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
MatrixMarket = "4d4711f2-db25-561a-b6b3-d35e7d4047d3"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
BenchmarkTools = "~1.2.0"
MatrixMarket = "~0.3.1"
PlutoUI = "~0.7.19"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0bc60e3006ad95b4bb7497698dd7c6d649b9bc06"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "61adeb0823084487000600ef8b1c00cc2474cd47"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.2.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MatrixMarket]]
deps = ["Compat", "SparseArrays"]
git-tree-sha1 = "54d39ccb57d29aefa666418bca8ca5598ebd8225"
uuid = "4d4711f2-db25-561a-b6b3-d35e7d4047d3"
version = "0.3.1"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "e071adf21e165ea0d904b595544a8e514c8bb42c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.19"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═6339ca2a-fb15-41fe-b3e3-5e7fbbe5095f
# ╠═a15f08d6-8b69-4ed0-b22b-f920f98db371
# ╠═e58d2a23-799c-44fa-9956-9ecc27e4fc8f
# ╠═a9a46efc-aabf-4cb9-8980-cc74b245ace5
# ╠═43b6eb91-e23c-4a2a-ab0a-5c219e9a965f
# ╠═0b87b6c3-b074-4f7f-8016-725425091e7e
# ╠═c8e428df-a6ce-48de-a64a-ad0b0c9077d3
# ╠═8b6a63e7-5f92-4659-8dcb-239e29dc65b1
# ╠═b6c1cb40-6bb2-4fca-ba94-dd77986c8b13
# ╠═dc12e43f-cb4c-465d-9da3-e4398ca041a4
# ╠═fa133740-2aca-4948-87cc-45f2941e608f
# ╠═89a843ef-4eb6-4472-ab7a-f9b85e3c7200
# ╠═f2574652-9f34-4fb0-aed0-b5e09807f087
# ╠═5d781032-288d-4bfc-826b-31afdb8040ea
# ╠═702db5d8-7d85-41e6-ae14-1b3fbb84f7c9
# ╠═a8bc2ff1-7ab0-4c6b-bf5b-31e43458ded0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
