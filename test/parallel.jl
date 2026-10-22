# The verification suite in several processes.
#
#     julia --startup-file=no --threads=8 test/parallel.jl [groups...]
#
# Most of the run time of the suite is compilation of test code, which one process
# cannot share between threads. This driver gives the groups to ALLSPICE_TEST_JOBS
# processes (default: a quarter of the CPU threads, 8 at most), each of which runs `runtests.jl` with a part of the groups and
# with the thread count of the driver. The driver prints the output of each process
# when it ends, then the time of each group and the wall time of the whole run, and it
# exits with status 1 if a process failed.
#
# The driver balances the processes with the group times of the last run, which it
# keeps in `test/.group_seconds` (not under version control). Without that file all
# groups have the same weight, so the second run is faster than the first.

const TESTDIR = @__DIR__
include(joinpath(TESTDIR, "kernels.jl"))
const TIMES_FILE = joinpath(TESTDIR, ".group_seconds")
const TIMES_HEADER = "Time of each group (s), longest first"

"Group names of `runtests.jl`, in the order of the sequential run."
function all_groups()
    text = read(joinpath(TESTDIR, "runtests.jl"), String)
    m = match(r"const GROUPS = \[(.*?)\]"s, text)
    m === nothing && error("no GROUPS list in runtests.jl")
    groups = [String(x.captures[1]) for x in eachmatch(r"\"([a-z0-9_]+)\"", m.captures[1])]
    @assert !isempty(groups) && allunique(groups)
    return groups
end

"Group times from `text` (the table that `runtests.jl` prints). Postcondition: all times >= 0."
function parse_times(text::AbstractString)
    times = Dict{String,Float64}()
    at = findfirst(TIMES_HEADER, text)
    at === nothing && return times
    for line in split(text[last(at):end], '\n')[2:end]
        m = match(r"^  ([a-z0-9_]+)\s+([0-9.]+)$", line)
        m === nothing && break
        m.captures[1] == "all" || (times[m.captures[1]] = parse(Float64, m.captures[2]))
    end
    @assert all(>=(0), values(times))
    return times
end

"""
    partition(groups, weight, njobs)

Give each group to one of `njobs` lists, the heaviest group first to the lightest list.
Each group can run alone, so any division is valid. Postcondition: the lists hold each
group once, each in the order of `groups`.
"""
function partition(groups::Vector{String}, weight, njobs::Integer)
    @assert njobs >= 1
    units = [[g] for g in groups]
    jobs = [String[] for _ in 1:min(njobs, length(units))]
    load = zeros(length(jobs))
    for u in sort(units; by=u -> sum(weight, u), rev=true)
        k = argmin(load)
        append!(jobs[k], u)
        load[k] += sum(weight, u)
    end
    order = Dict(g => i for (i, g) in enumerate(groups))
    foreach(j -> sort!(j; by=g -> order[g]), jobs)
    @assert sort(reduce(vcat, jobs)) == sort(groups)
    return jobs
end

function main(args)
    groups = all_groups()
    selected = isempty(args) ? groups : args
    unknown = setdiff(selected, groups)
    isempty(unknown) || error("unknown test group(s) $(unknown); the groups are $(groups)")
    selected = filter(in(selected), groups)
    njobs = parse(Int, get(ENV, "ALLSPICE_TEST_JOBS", string(clamp(Sys.CPU_THREADS ÷ 4, 1, 8))))
    known = isfile(TIMES_FILE) ? parse_times(TIMES_HEADER * "\n" * read(TIMES_FILE, String)) :
        Dict{String,Float64}()
    mean = isempty(known) ? 1.0 : sum(values(known)) / length(known)
    jobs = partition(selected, g -> get(known, g, mean), njobs)

    Kernels.fetch_kernels()     # once here, so that the processes do not download in parallel
    logdir = mktempdir()
    julia = `$(Base.julia_cmd()) --startup-file=no --threads=$(Threads.nthreads())
        --project=$TESTDIR $(joinpath(TESTDIR, "runtests.jl"))`
    start = time()
    procs = map(enumerate(jobs)) do (k, job)
        log = joinpath(logdir, "job$k.log")
        # One stream for both outputs: two opens of the path would overwrite each other.
        io = open(log, "w")
        (log, io) => run(pipeline(ignorestatus(`$julia $job`); stdout=io, stderr=io); wait=false)
    end
    times = Dict{String,Float64}()
    failed = Int[]
    for (k, ((log, io), proc)) in enumerate(procs)
        wait(proc)
        close(io)
        text = read(log, String)
        println("\n===== process ", k, " of ", length(procs), ": ", join(jobs[k], " "), " =====")
        print(text)
        merge!(times, parse_times(text))
        success(proc) || push!(failed, k)
    end
    wall = time() - start

    println("\n", TIMES_HEADER, ", all processes:")
    table = sort(collect(times); by=last, rev=true)
    lines = ["  " * rpad(g, 18) * lpad(round(t; digits=1), 7) for (g, t) in table]
    foreach(println, lines)
    println("  ", rpad("sum", 18), lpad(round(sum(values(times); init=0.0); digits=1), 7))
    println("\nWall time ", round(wall; digits=1), " s in ", length(procs), " processes")
    isempty(times) || write(TIMES_FILE, join(lines, "\n") * "\n")
    if !isempty(failed)
        printstyled("FAILED in process(es) ", failed, "\n"; color=:red, bold=true)
        exit(1)
    end
    return nothing
end

main(ARGS)
