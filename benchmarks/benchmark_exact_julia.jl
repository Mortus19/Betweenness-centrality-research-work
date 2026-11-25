# Exact shortest temporal betweenness benchmark for .patg graphs.
# Julia >= 1.10; install TemporalGraphs.jl once:
#   julia -e 'using Pkg; Pkg.add(url="https://github.com/aurorarossi/TemporalGraphs.jl")'
# Usage:
#   julia --startup-file=no --threads=4 benchmark_exact_julia.jl INPUT.patg VALUES.csv [SUMMARY.csv]
#
# Match the user's C++ preprocessing:
#   - skip the .patg header (one first line)
#   - discard self-loops and identical (u,v,t) edges
#   - renumber active original vertex IDs in ascending order
#   - use directed edges and strictly increasing timestamps
#   - use minimum hops, with no waiting-time restriction
#   - normalize by n*(n-1), including unreachable ordered pairs

using TemporalGraphs
using Printf

function read_cpp_compatible_patg(filename::AbstractString)
    raw_edges = NTuple{3, Int}[]
    open(filename, "r") do file
        eof(file) && error("Empty .patg file: $filename")
        readline(file)  # The C++ reader also discards the first line.
        for (lineno, line) in enumerate(eachline(file))
            line = strip(line)
            (isempty(line) || startswith(line, "#")) && continue
            parts = split(line)
            length(parts) >= 3 || error("Malformed edge at line $(lineno + 1): $line")
            u = parse(Int, parts[1])
            v = parse(Int, parts[2])
            t = parse(Int, parts[3])
            u == v && continue
            push!(raw_edges, (u, v, t))
        end
    end

    isempty(raw_edges) && error("No non-loop temporal edges in $filename")
    sort!(raw_edges)
    unique!(raw_edges)

    original_ids = Int[]
    sizehint!(original_ids, 2 * length(raw_edges))
    for (u, v, _) in raw_edges
        push!(original_ids, u)
        push!(original_ids, v)
    end
    sort!(original_ids)
    unique!(original_ids)
    n = length(original_ids)
    n >= 2 || error("At least two active vertices are required")

    id_to_index = Dict(id => i for (i, id) in enumerate(original_ids))
    # Transition time 1 means edges (u,v,t) and (v,w,t') can follow
    # each other exactly when t' > t (timestamps are integers).
    # +1 reproduces the harmless timestamp shift in the C++ reader.
    converted_edges = Vector{NTuple{4, Int}}(undef, length(raw_edges))
    for (i, (u, v, t)) in enumerate(raw_edges)
        converted_edges[i] = (id_to_index[u], id_to_index[v], Base.checked_add(t, 1), 1)
    end

    graph = OrderedEdgeList(n, converted_edges)
    return graph, original_ids, length(raw_edges)
end

function exact_shortest(g)
    return temporal_betweenness(g, MinimumHops())
end

function csv_field(value)
    s = string(value)
    return "\"" * replace(s, "\"" => "\"\"") * "\""
end

function benchmark(input_path::String, values_path::String, summary_path::Union{Nothing, String})
    isfile(input_path) || error("Graph not found: $input_path")

    # Preprocessing is measured independently of the shortest-BC algorithm.
    t0 = time_ns()
    graph, original_ids, m = read_cpp_compatible_patg(input_path)
    preprocess_s = (time_ns() - t0) / 1e9
    n = length(original_ids)

    # Warm up the SAME algorithm on a tiny graph to trigger Julia JIT compilation.
    # This is deliberately not included in compute_s.
    warmup_graph = OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 2, 1)])
    exact_shortest(warmup_graph)
    GC.gc()

    t1 = time_ns()
    bc_raw = exact_shortest(graph)
    compute_s = (time_ns() - t1) / 1e9
    length(bc_raw) == n || error("Wrong result length: $(length(bc_raw)) instead of $n")

    factor = 1.0 / (Float64(n) * Float64(n - 1))
    bc_normalized = Float64.(bc_raw) .* factor

    mkpath(dirname(abspath(values_path)))
    open(values_path, "w") do io
        println(io, "index,original_id,bc_raw,bc_normalized")
        for i in 1:n
            @printf(io, "%d,%d,%.17g,%.17g\n", i - 1, original_ids[i], Float64(bc_raw[i]), bc_normalized[i])
        end
    end

    @printf("Dataset           : %s\n", basename(input_path))
    @printf("Algorithm         : Julia / TemporalGraphs.jl / MinimumHops (strict)\n")
    @printf("Vertices          : %d\n", n)
    @printf("Temporal edges    : %d\n", m)
    @printf("Julia threads     : %d\n", Threads.nthreads())
    @printf("Preprocessing (s) : %.6f\n", preprocess_s)
    @printf("Compute time (s)  : %.6f\n", compute_s)
    @printf("Total (s)         : %.6f (excluding Julia startup and JIT warmup)\n", preprocess_s + compute_s)
    @printf("BC sum normalized : %.12g\n", sum(bc_normalized))
    println("Values saved     : ", abspath(values_path))

    if summary_path !== nothing
        mkpath(dirname(abspath(summary_path)))
        write_header = !isfile(summary_path) || filesize(summary_path) == 0
        open(summary_path, "a") do io
            if write_header
                println(io, "dataset,algorithm,julia_version,threads,n,m,preprocess_s,compute_s,total_s,values_file")
            end
            @printf(io, "%s,%s,%s,%d,%d,%d,%.9f,%.9f,%.9f,%s\n",
                    csv_field(basename(input_path)),
                    csv_field("TemporalGraphs.jl MinimumHops strict"),
                    csv_field(VERSION),
                    Threads.nthreads(), n, m,
                    preprocess_s, compute_s, preprocess_s + compute_s,
                    csv_field(abspath(values_path)))
        end
        println("Summary appended : ", abspath(summary_path))
    end
end

function main()
    if length(ARGS) < 2 || length(ARGS) > 3
        println(stderr, "Usage: julia --threads=4 benchmark_exact_julia.jl INPUT.patg VALUES.csv [SUMMARY.csv]")
        exit(2)
    end
    summary = length(ARGS) == 3 ? ARGS[3] : nothing
    benchmark(ARGS[1], ARGS[2], summary)
end

main()
