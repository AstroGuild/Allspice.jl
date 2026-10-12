# Writing an EK file, part 1: pages and trees.
#
# This file ports the page allocator of CSPICE (`zzekpage.c`) and the insert and delete
# operations of the EK B*-tree (`zzektr*.c`). The layouts are those that `file.jl`
# describes. A write handle holds the three DAS address spaces in memory, so a "page
# read" of CSPICE is a copy of 256 integers and a "page write" copies them back.
#
# The port follows the C source statement by statement, because the words that a file
# holds after a sequence of calls must be the words that CSPICE writes: the order of
# the page allocations, the splits and the rotations, and also the words that CSPICE
# leaves in the unused part of a page. For that reason a tree page is handled through
# `_ekw_Page`, whose indexes are the 0-based offsets of the f2c source (`page[4]` is the
# fifth integer of the page), and a statement here can be compared with its C statement
# one to one.
#
# Nothing here keeps state outside the handle.

"Data words of a page of type 1 (character), 2 (double), 3 (integer)."
const _ekw_DSIZE = (_EK_CPSIZE, _EK_DPSIZE, _EK_IPSIZE)

# Positions of the paging metadata in the first 256 integers of the file (`ekfilpar.inc`):
# for each type the number of pages, the number of free pages, the head of the free list.
const _ekw_NPAGE = (4, 9, 14)
const _ekw_NFREE = (5, 10, 15)
const _ekw_FREEP = (6, 11, 16)

"""
    _ekw_check(h, who) -> Nothing

The check of `zzekpgch` for write access: `h` is open (`:dasnosuchhandle`), writable
(`:dasinvalidaccess`) and holds a paged EK (`:wrongarchitecture`, `:invalidformat`).
"""
@inline function _ekw_check(h::DASHandle, who::AbstractString)
    _das_check_write(h, who)
    _ek_check(h, h.path, who)
    return nothing
end

"Address before the first word of page `p` of type `t` (`zzekpgbs`)."
@inline _ekw_base(t::Int, p::Integer) =
    t == 1 ? _ek_cbase(p) : t == 2 ? _ek_dbase(p) : _ek_ibase(p)

"Page of type `t` that holds the address `addr`, and its base (`zzekpgpg`)."
@inline function _ekw_pgpg(t::Int, addr::Integer)
    a = Int(addr)
    p = t == 1 ? (a + 1023) ÷ 1024 : t == 2 ? (a + 127) ÷ 128 : (a - 1) ÷ 256
    return p, _ekw_base(t, p)
end

"""
    _ekw_sei(h, addr, val) -> Nothing

Write `val` as five base-128 characters at the character address `addr` (`zzeksei`,
`prtenc`). Precondition: `0 <= val < 128^5`, else `SpiceError(:outofrange)`.
"""
function _ekw_sei(h::DASHandle, addr::Integer, val::Integer)
    (0 <= val < 34359738368) || throw(SpiceError(:outofrange,
        "the integer $val cannot be encoded in five base-128 characters"))
    n = Int(val)
    chars = h.chars
    for k in 4:-1:0                      # invariant: digits after k are stored
        chars[addr+k] = (n % 128) % UInt8
        n ÷= 128
    end
    return nothing
end

"Write the bytes `src[first:first+n-1]` at the character address `addr` (`dasudc`)."
@inline function _ekw_udc(h::DASHandle, addr::Integer, src::AbstractVector{UInt8}, first::Int,
        n::Int)
    n > 0 && copyto!(h.chars, Int(addr), src, first, n)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Page allocation
# ---------------------------------------------------------------------------------------

"""
    _ekw_pgin(h) -> Nothing

Set up the paging metadata of a new, empty DAS file (`zzekpgin`). Precondition: the
file holds no data, else `SpiceError(:dasnotempty)`. Postcondition: the first 256
integers hold the architecture code, the page sizes and the page counts (all zero).
"""
function _ekw_pgin(h::DASHandle)
    _das_check_write(h, "zzekpgin")
    (isempty(h.chars) && isempty(h.dps) && isempty(h.ints)) || throw(SpiceError(:dasnotempty,
        "File $(h.path) contains data; LASTC = $(length(h.chars)); LASTD = " *
        "$(length(h.dps)); LASTI = $(length(h.ints))."))
    append!(h.ints, zeros(Int32, 256))
    ints = h.ints
    ints[1] = _EK_ARCH
    ints[2] = 1024
    ints[7] = 128
    ints[12] = 256
    ints[3] = 0
    ints[8] = 0
    ints[13] = 256
    return nothing
end

"""
    _ekw_pgan(h, t) -> (p, base)

Append a new page of type `t` to the file (`zzekpgan`): 1024 blanks, 128 zeros or 256
zeros. Postcondition: the page count of the type is one more, and `p` is that count.
"""
function _ekw_pgan(h::DASHandle, t::Int)
    ints = h.ints
    if t == 1
        n = length(h.chars)
        resize!(h.chars, n + 1024)
        fill!(view(h.chars, n+1:n+1024), UInt8(' '))
    elseif t == 2
        n = length(h.dps)
        resize!(h.dps, n + 128)
        fill!(view(h.dps, n+1:n+128), 0.0)
    elseif t == 3
        n = length(ints)
        resize!(ints, n + 256)
        fill!(view(ints, n+1:n+256), Int32(0))
    else
        throw(SpiceError(:invalidtype, "The data type code $t was not recognized."))
    end
    np = Int(ints[_ekw_NPAGE[t]])
    ints[_ekw_NPAGE[t]] = np + 1
    return np + 1, _ekw_base(t, np + 1)
end

"""
    _ekw_pgal(h, t) -> (p, base)

Allocate a page of type `t` (`zzekpgal`): the head of the free list of the type if the
list has a page, else a new page. A page from the free list keeps its old words, as in
CSPICE; its first word held the link of the free list.
"""
function _ekw_pgal(h::DASHandle, t::Int)
    (1 <= t <= 3) || throw(SpiceError(:invalidtype, "The data type code $t was not recognized."))
    ints = h.ints
    p = Int(ints[_ekw_FREEP[t]])
    p > 0 || return _ekw_pgan(h, t)
    base = _ekw_base(t, p)
    next = t == 1 ? _ek_gei(h.chars, base + 1) : t == 2 ? round(Int, h.dps[base+1]) :
           Int(ints[base+1])
    ints[_ekw_NFREE[t]] -= 1
    ints[_ekw_FREEP[t]] = next
    return p, base
end

"""
    _ekw_pgfr(h, t, p) -> Nothing

Put page `p` of type `t` at the head of the free list of its type (`zzekpgfr`).
Precondition: the page exists, else `SpiceError(:invalidindex)`. Only the first word
of the page changes: it takes the old head of the list.
"""
function _ekw_pgfr(h::DASHandle, t::Int, p::Integer)
    (1 <= t <= 3) || throw(SpiceError(:invalidtype, "The data type code $t was not recognized."))
    ints = h.ints
    np = Int(ints[_ekw_NPAGE[t]])
    (1 <= p <= np) || throw(SpiceError(:invalidindex, "Attempt to free non-existent " *
        "$(("CHR", "DP", "INT")[t]) page. Page number = $p; valid range is 1:$np"))
    head = Int(ints[_ekw_FREEP[t]])
    base = _ekw_base(t, p)
    if t == 1
        _ekw_sei(h, base + 1, head)
    elseif t == 2
        h.dps[base+1] = Float64(head)
    else
        ints[base+1] = head
    end
    ints[_ekw_FREEP[t]] = p
    ints[_ekw_NFREE[t]] += 1
    return nothing
end

"Link count of page `p` of type `t` (`zzekglnk`)."
@inline function _ekw_glnk(h::DASHandle, t::Int, p::Integer)
    base = _ekw_base(t, p)
    return t == 1 ? _ek_gei(h.chars, base + 1020) : t == 2 ? round(Int, h.dps[base+128]) :
           Int(h.ints[base+256])
end

"Set the link count of page `p` of type `t` to `n` (`zzekslnk`)."
@inline function _ekw_slnk(h::DASHandle, t::Int, p::Integer, n::Integer)
    base = _ekw_base(t, p)
    if t == 1
        _ekw_sei(h, base + 1020, n)
    elseif t == 2
        h.dps[base+128] = Float64(n)
    else
        h.ints[base+256] = n
    end
    return nothing
end

"Forward pointer of page `p` of type `t` (`zzekgfwd`)."
@inline function _ekw_gfwd(h::DASHandle, t::Int, p::Integer)
    base = _ekw_base(t, p)
    return t == 1 ? _ek_gei(h.chars, base + 1015) : t == 2 ? round(Int, h.dps[base+127]) :
           Int(h.ints[base+255])
end

"Set the forward pointer of page `p` of type `t` to the page `fwd` (`zzeksfwd`)."
@inline function _ekw_sfwd(h::DASHandle, t::Int, p::Integer, fwd::Integer)
    base = _ekw_base(t, p)
    if t == 1
        _ekw_sei(h, base + 1015, fwd)
    elseif t == 2
        h.dps[base+127] = Float64(fwd)
    else
        h.ints[base+255] = fwd
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Integer page buffers with the 0-based offsets of the C source
# ---------------------------------------------------------------------------------------

"""
    _ekw_Page

A copy of one integer page (256 integers). `page[k]` is the integer at the 0-based
offset `k` of the page, as in the f2c source of CSPICE; values are read as `Int`.
"""
struct _ekw_Page
    v::Vector{Int32}
end
_ekw_Page() = _ekw_Page(zeros(Int32, 256))
@inline Base.getindex(p::_ekw_Page, k::Integer) = Int(p.v[k+1])
@inline Base.setindex!(p::_ekw_Page, x::Integer, k::Integer) = (p.v[k+1] = x)

"Throws `SpiceError(:invalidindex)` unless `p` is an integer page of the file."
@inline function _ekw_ipage_check(h::DASHandle, p::Integer)
    np = Int(h.ints[14])
    (1 <= p <= np) || throw(SpiceError(:invalidindex, "INT page = $p; valid range is [1:$np]"))
    return nothing
end

"Read integer page `p` into the buffer `page` (`zzekpgri`); returns `page`."
function _ekw_pgri!(page::_ekw_Page, h::DASHandle, p::Integer)
    _ekw_ipage_check(h, p)
    copyto!(page.v, 1, h.ints, _ek_ibase(p) + 1, 256)
    return page
end

"A copy of integer page `p` (`zzekpgri`)."
_ekw_pgri(h::DASHandle, p::Integer) = _ekw_pgri!(_ekw_Page(Vector{Int32}(undef, 256)), h, p)

"Write the buffer `page` to integer page `p` (`zzekpgwi`)."
function _ekw_pgwi(h::DASHandle, p::Integer, page::_ekw_Page)
    _ekw_ipage_check(h, p)
    copyto!(h.ints, _ek_ibase(p) + 1, page.v, 1, 256)
    return nothing
end

"Copy `n` integers from offset `so` of `src` to offset `to` of `dst` (`movei`); `src !== dst`."
@inline function _ekw_move!(dst::_ekw_Page, to::Int, src::_ekw_Page, so::Int, n::Int)
    n > 0 && copyto!(dst.v, to + 1, src.v, so + 1, n)
    return nothing
end

"Set `n` integers from offset `k` of `page` to zero (`cleari`); nothing if `n <= 0`."
@inline function _ekw_clear!(page::_ekw_Page, k::Int, n::Int)
    n > 0 && fill!(view(page.v, k+1:k+n), Int32(0))
    return nothing
end

"Integer at the 0-based offset `k` of the page with the base `b`, read from the file."
@inline _ekw_w(ints::Vector{Int32}, b::Int, k::Int) = Int(ints[b+k+1])

# ---------------------------------------------------------------------------------------
# Tree: lookups that the writers need
# ---------------------------------------------------------------------------------------

"Number of keys of the node `node` of the tree `tree` (`zzektrnk`)."
@inline _ekw_trnk(h::DASHandle, tree::Integer, node::Integer) =
    Int(h.ints[_ek_ibase(node)+(tree == node ? 5 : 1)])

"""
    _ekw_trpi(h, tree, key) -> (parent, pkey, poffst, lpidx, lpkey, lsib, rpidx, rpkey, rsib)

Parent information of the node that holds `key` (`zzektrpi`): the parent page, the first
key of the parent (absolute), the key offset of the parent, and for the keys of the
parent that are left and right of the node their position, their value in the parent and
the sibling of the node on that side. All are 0 if the key is in the root; the values of
a side are 0 if the node is the first or the last child. Precondition: `1 <= key <=`
number of keys, else `SpiceError(:indexoutofrange)`; `:itemnotfound` for a damaged tree.
"""
function _ekw_trpi(h::DASHandle, tree::Integer, key::Integer)
    ints = h.ints
    lkey = Int(key)
    b = _ek_ibase(tree)
    parent = pkey = poffst = lpidx = lpkey = lsib = rpidx = rpkey = rsib = 0
    totkey = _ekw_w(ints, b, 2)
    (1 <= lkey <= totkey) || throw(SpiceError(:indexoutofrange,
        "Key = $lkey; valid range = 1:$totkey. Tree = $tree, file = $(h.path)"))
    prev = _ek_lstlei(ints, b + 6, _ekw_w(ints, b, 4), lkey)
    prvkey = prev > 0 ? _ekw_w(ints, b, prev + 4) : 0
    prvkey == lkey && return (0, 0, 0, 0, 0, 0, 0, 0, 0)
    offset = prvkey
    parent = Int(tree)
    pkey = _ekw_w(ints, b, 5)
    maxkey = _ekw_w(ints, b, 4)
    if prev > 0
        lpidx = prev
        lpkey = _ekw_w(ints, b, lpidx + 4)
        lsib = _ekw_w(ints, b, lpidx + 87)
    end
    if prev < maxkey
        rpidx = prev + 1
        rpkey = _ekw_w(ints, b, rpidx + 4)
        rsib = _ekw_w(ints, b, rpidx + 88)
    end
    child = _ekw_w(ints, b, prev + 88)
    found = false
    hops = 0
    while child > 0 && !found
        # invariant: the key is in the subtree of `child`, whose keys start after `offset`;
        # the nine results describe the node above `child`
        hops += 1
        hops > _EK_TRMXDP && break
        b = _ek_ibase(child)
        newkey = lkey - offset
        prev = _ek_lstlei(ints, b + 2, _ekw_w(ints, b, 0), newkey)
        prvkey = prev > 0 ? _ekw_w(ints, b, prev) : 0
        if prvkey == newkey
            found = true
        else
            parent = child
            poffst = offset
            pkey = _ekw_w(ints, b, 1) + offset
            maxkey = _ekw_w(ints, b, 0)
            if prev > 0
                lpidx = prev
                lpkey = _ekw_w(ints, b, lpidx)
                lsib = _ekw_w(ints, b, lpidx + 63)
            else
                lpidx = lpkey = lsib = 0
            end
            if prev < maxkey
                rpidx = prev + 1
                rpkey = _ekw_w(ints, b, rpidx)
                rsib = _ekw_w(ints, b, rpidx + 64)
            else
                rpidx = rpkey = rsib = 0
            end
            offset = prvkey + offset
            child = _ekw_w(ints, b, prev + 64)
        end
    end
    found || throw(SpiceError(:itemnotfound, "Key $lkey; valid range = 1:$totkey. Tree = " *
        "$tree, file = $(h.path).  Key was not found.  This probably indicates a corrupted " *
        "file or a bug in the EK code."))
    return (parent, pkey, poffst, lpidx, lpkey, lsib, rpidx, rpkey, rsib)
end

"""
    _ekw_trki(h, tree, nodkey, n) -> key

The absolute value of key number `n` of the node that holds the key `nodkey`
(`zzektrki`). Precondition: `0 <= n <=` number of keys of the node (`:invalidindex`).
"""
function _ekw_trki(h::DASHandle, tree::Integer, nodkey::Integer, n::Integer)
    _, node, noffst, _, level = _ek_trlk(h.ints, tree, nodkey)
    size = _ekw_trnk(h, tree, node)
    (0 <= n <= size) || throw(SpiceError(:invalidindex,
        "Key index = $n; valid range in node $node is 1:$size"))
    base = _ek_ibase(node)
    return Int(h.ints[base+(level == 1 ? 5 : 1)+n]) + noffst
end

"""
    _ekw_trsb(h, tree, key) -> (lsib, lkey, rsib, rkey)

The left and the right sibling of the node that holds `key`, and the first key
(absolute) of each (`zzektrsb`); 0 where the node has no sibling on that side.
"""
function _ekw_trsb(h::DASHandle, tree::Integer, key::Integer)
    ints = h.ints
    parent, _, poffst, lpidx, _, lsib, rpidx, rpkey, rsib = _ekw_trpi(h, tree, key)
    keybas = parent == tree ? 5 : 1
    lkey = 0
    rkey = 0
    if lpidx > 1
        loffst = Int(ints[_ek_ibase(parent)+keybas+lpidx-1]) + poffst
        lkey = Int(ints[_ek_ibase(lsib)+2]) + loffst
    elseif lpidx == 1
        lkey = Int(ints[_ek_ibase(lsib)+2]) + poffst
    end
    if rpidx > 0
        rkey = Int(ints[_ek_ibase(rsib)+2]) + rpkey + poffst
    end
    return lsib, lkey, rsib, rkey
end

"Key of the first entry of the tree `tree` whose value is `ival`, or 0 (`zzektrls`)."
function _ekw_trls(h::DASHandle, tree::Integer, ival::Integer)
    vals = _ek_tree_values(h.ints, tree)
    for i in eachindex(vals)             # invariant: the entries before i have another value
        vals[i] == ival && return i
    end
    return 0
end

@noinline _ekw_bug(msg::AbstractString) = throw(SpiceError(:bug, msg))

# ---------------------------------------------------------------------------------------
# Tree: creation and insertion
# ---------------------------------------------------------------------------------------

"""
    _ekw_trit(h) -> tree

Make an empty tree and return its root page (`zzektrit`). Postcondition: the root page
is zero but for the tree type 1, the node count 1 and the depth 1.
"""
function _ekw_trit(h::DASHandle)
    p, base = _ekw_pgal(h, 3)
    ints = h.ints
    fill!(view(ints, base+1:base+256), Int32(0))
    ints[base+_EK_TRTYPE] = 1
    ints[base+_EK_TRNNOD] = 1
    ints[base+_EK_TRDPTH] = 1
    return p
end

"""
    _ekw_trui(h, tree, key, value) -> overflow

Insert `value` at the key `key` without balancing (`zzektrui`): the keys from `key` on
move up by one. The result tells that the node that took the key now holds one key more
than its capacity. Precondition: `1 <= key <= nkeys + 1`, else `:invalidindex`.
"""
function _ekw_trui(h::DASHandle, tree::Integer, key::Int, value::Integer)
    root = Int(tree)
    rpage = _ekw_pgri(h, root)
    totkey = rpage[2]
    (1 <= key <= totkey + 1) || throw(SpiceError(:invalidindex,
        "Key = $key. Valid range is 1:$(totkey + 1).  File = $(h.path)."))
    nnode = rpage[1]
    depth = rpage[3]
    if nnode == 1
        nkeys = totkey
        rpage[2] = nkeys + 1
        rpage[4] = nkeys + 1
        for i in nkeys:-1:key
            rpage[i+5] = rpage[i+4] + 1
            rpage[i+172] = rpage[i+171]
        end
        for i in nkeys+1:-1:key
            rpage[i+88] = rpage[i+87]
        end
        rpage[key+4] = key
        rpage[key+171] = value
        rpage[key+87] = 0
        nkeys += 1
        _ekw_pgwi(h, root, rpage)
        return nkeys == 83
    elseif key == totkey + 1
        _, target, _, _, _ = _ek_trlk(h.ints, tree, key - 1)
        tpage = _ekw_pgri(h, target)
        nkeys = tpage[0]
        keyloc = nkeys + 2
        datloc = nkeys + 129
        kidloc = nkeys + 65
        tpage[keyloc-1] = tpage[keyloc-2] + 1
        tpage[datloc-1] = value
        tpage[kidloc-1] = 0
        tpage[0] = tpage[0] + 1
        rpage[2] = totkey + 1
        _ekw_pgwi(h, root, rpage)
        _ekw_pgwi(h, target, tpage)
        return nkeys == 62
    end
    _, target, _, next, level = _ek_trlk(h.ints, tree, key)
    local tpage::_ekw_Page, tnkeys::Int, parent::Int, pkey::Int, rpidx::Int
    if level == depth
        parent, pkey, _, _, _, _, rpidx, _, _ = _ekw_trpi(h, tree, key)
        tpage = _ekw_pgri(h, target)
        tnkeys = tpage[0]
        tnkeys > 62 && throw(SpiceError(:nodetoofull, "Node = $target. Tree = $tree. File = " *
            "$(h.path). Key count = $tnkeys; max allowed, including overflow, is 63."))
        for i in tnkeys:-1:next
            tpage[i+1] = tpage[i] + 1
        end
        for i in tnkeys:-1:next
            tpage[i+128] = tpage[i+127]
        end
        for i in tnkeys+1:-1:next
            tpage[i+64] = tpage[i+63]
        end
        tpage[next+127] = value
    else
        # The key is not in a leaf: its predecessor is the last key of a leaf, and the
        # new key goes after that one.
        _, target, _, prev, level = _ek_trlk(h.ints, tree, key - 1)
        parent, pkey, _, _, _, _, rpidx, _, _ = _ekw_trpi(h, tree, key - 1)
        tpage = _ekw_pgri(h, target)
        tnkeys = tpage[0]
        tnkeys > 63 && throw(SpiceError(:nodetoofull, "Node = $target. Tree = $tree. File = " *
            "$(h.path). Key count = $tnkeys; max allowed, including overflow, is 63."))
        tpage[prev+1] = prev + 1
        tpage[prev+128] = value
        tpage[prev+65] = 0
    end
    tpage[0] = tnkeys + 1
    overfl = tnkeys == 62
    _ekw_pgwi(h, target, tpage)
    while parent != root
        # invariant: the nodes below `parent` on the path hold the new key; in `parent`
        # the keys from position `rpidx` on are those after the new key
        paren2, pkey2, _, _, _, _, rpidx2, _, _ = _ekw_trpi(h, tree, pkey)
        if rpidx > 0
            _ekw_pgri!(tpage, h, parent)
            for i in rpidx:tpage[0]
                tpage[i] = tpage[i] + 1
            end
            _ekw_pgwi(h, parent, tpage)
        end
        parent = paren2
        pkey = pkey2
        rpidx = rpidx2
    end
    if rpidx > 0
        for i in rpidx:rpage[4]
            rpage[i+4] = rpage[i+4] + 1
        end
    end
    rpage[2] = totkey + 1
    _ekw_pgwi(h, root, rpage)
    return overfl
end

"""
    _ekw_tr13(h, tree) -> Nothing

Split the root that holds 83 keys into a root with one key and two children with 41
keys each (`zzektr13`). The depth of the tree grows by one.
"""
function _ekw_tr13(h::DASHandle, tree::Integer)
    root = Int(tree)
    rpage = _ekw_pgri(h, root)
    nrkeys = rpage[4]
    nrkeys == 83 || _ekw_bug("Number of keys in root = $nrkeys; should be 83.")
    child1, _ = _ekw_pgal(h, 3)
    child2, _ = _ekw_pgal(h, 3)
    c1page = _ekw_Page()
    c1page[0] = 41
    _ekw_move!(c1page, 1, rpage, 5, 41)
    _ekw_move!(c1page, 128, rpage, 172, 41)
    _ekw_move!(c1page, 64, rpage, 88, 42)
    c2page = _ekw_Page()
    c2page[0] = 41
    middle = 42
    _ekw_move!(c2page, 1, rpage, middle + 5, 41)
    _ekw_move!(c2page, 128, rpage, middle + 172, 41)
    _ekw_move!(c2page, 64, rpage, middle + 88, 42)
    delta = rpage[middle+4]
    for i in 1:41
        c2page[i] = c2page[i] - delta
    end
    rpage[5] = rpage[middle+4]
    rpage[172] = rpage[middle+171]
    rpage[88] = child1
    rpage[89] = child2
    rpage[4] = 1
    rpage[3] = rpage[3] + 1
    rpage[1] = rpage[1] + 2
    _ekw_clear!(rpage, 6, 82)
    _ekw_clear!(rpage, 173, 82)
    _ekw_clear!(rpage, 90, 82)
    _ekw_pgwi(h, root, rpage)
    _ekw_pgwi(h, child1, c1page)
    _ekw_pgwi(h, child2, c2page)
    return nothing
end

"Key, data and child bases (0-based offsets) of a node: the root or another node."
@inline _ekw_bases(isroot::Bool) = isroot ? (5, 172, 88) : (1, 128, 64)

"Throws `:bug` unless the children `pkidx` and `pkidx + 1` of `ppage` are `left` and `right`."
@inline function _ekw_sibling_check(ppage::_ekw_Page, kidbas::Int, pkidx::Int, left::Int,
        right::Int, parent::Int)
    lsib = ppage[kidbas+pkidx-1]
    rsib = ppage[kidbas+pkidx]
    (lsib == left && rsib == right) || _ekw_bug("LEFT, RIGHT, PARENT, and PKIDX are " *
        "inconsistent. LEFT = $left; RIGHT = $right; PARENT = $parent; PKIDX = $pkidx; " *
        "LSIB derived from PARENT = $lsib; RSIB = $rsib.")
    return nothing
end

"Add `delta` to the node count of the tree, in `ppage` if the parent is the root (`zzektr23`, `zzektr32`)."
@inline function _ekw_nnode!(h::DASHandle, root::Int, parent::Int, ppage::_ekw_Page, delta::Int)
    if parent == root
        ppage[1] = ppage[1] + delta
    else
        h.ints[_ek_ibase(root)+2] += delta
    end
    return nothing
end

"""
    _ekw_tr23(h, tree, left, right, parent, pkidx) -> overflow

Split the two siblings `left` and `right`, which hold 125 keys together, into three
nodes with 41, 41 and 42 keys (`zzektr23`); the new node goes between them. `pkidx` is
the position in `parent` of the key between the siblings. The result tells that the
parent now holds one key more than its capacity.
"""
function _ekw_tr23(h::DASHandle, tree::Integer, left::Int, right::Int, parent::Int, pkidx::Int)
    lsize = 41
    msize = 41
    rsize = 42
    c1page = _ekw_pgri(h, left)
    c2page = _ekw_pgri(h, right)
    ppage = _ekw_pgri(h, parent)
    root = Int(tree)
    keybas, datbas, kidbas = _ekw_bases(parent == root)
    _ekw_sibling_check(ppage, kidbas, pkidx, left, right, parent)
    nlkeys = c1page[0]
    nrkeys = c2page[0]
    nlkeys + nrkeys == 125 || _ekw_bug("Number of keys in LEFT = $nlkeys; number of keys " *
        "in right = $nrkeys; but sum should be 125.")
    new, _ = _ekw_pgal(h, 3)
    c3page = _ekw_Page()
    lshift = nlkeys - (lsize + 1)
    ldelta = c1page[lsize+1]
    for i in 1:lshift
        c3page[i] = c1page[lsize+2+i-1] - ldelta
    end
    _ekw_move!(c3page, 128, c1page, lsize + 129, lshift)
    _ekw_move!(c3page, 64, c1page, lsize + 65, lshift + 1)
    ltrsiz = pkidx == 1 ? ppage[keybas] - 1 : ppage[keybas+pkidx-1] - ppage[keybas+pkidx-2] - 1
    lnmove = ltrsiz - ldelta
    c3page[lshift+1] = lnmove + 1
    c3page[lshift+128] = ppage[datbas+pkidx-1]
    rshift = nrkeys - (rsize + 1)
    rdelta = lnmove + 1
    for i in 1:rshift
        c3page[lshift+2+i-1] = c2page[i] + rdelta
    end
    _ekw_move!(c3page, lshift + 129, c2page, 128, rshift)
    _ekw_move!(c3page, lshift + 65, c2page, 64, rshift + 1)
    rnmove = c2page[rshift+1] - 1
    c3page[0] = msize
    npkeys = parent == root ? ppage[4] : ppage[0]
    for i in npkeys:-1:pkidx+1
        ppage[keybas+i] = ppage[keybas+i-1]
    end
    for i in npkeys:-1:pkidx+1
        ppage[datbas+i] = ppage[datbas+i-1]
    end
    for i in npkeys+1:-1:pkidx+1
        ppage[kidbas+i] = ppage[kidbas+i-1]
    end
    ppage[datbas+pkidx-1] = c1page[lsize+128]
    ppage[keybas+pkidx-1] = pkidx == 1 ? ldelta : ppage[keybas+pkidx-2] + ldelta
    ppage[datbas+pkidx] = c2page[rshift+128]
    ppage[keybas+pkidx] = ppage[keybas+pkidx-1] + lnmove + rnmove + 2
    ppage[kidbas+pkidx] = new
    local overfl::Bool
    if parent == root
        ppage[4] = ppage[4] + 1
        overfl = ppage[4] == 83
    else
        ppage[0] = ppage[0] + 1
        overfl = ppage[0] == 63
    end
    c1page[0] = lsize
    _ekw_clear!(c1page, lsize + 1, 63 - lsize)
    _ekw_clear!(c1page, lsize + 128, 63 - lsize)
    _ekw_clear!(c1page, lsize + 65, 64 - (lsize + 1))
    for i in 1:rsize
        c2page[i] = c2page[rshift+2+i-1] - (rnmove + 1)
    end
    for i in 1:rsize
        c2page[i+127] = c2page[rshift+129+i-1]
    end
    for i in 1:rsize+1
        c2page[i+63] = c2page[rshift+65+i-1]
    end
    c2page[0] = rsize
    _ekw_clear!(c2page, rsize + 1, 63 - rsize)
    _ekw_clear!(c2page, rsize + 128, 63 - rsize)
    _ekw_clear!(c2page, rsize + 65, 64 - (rsize + 1))
    _ekw_nnode!(h, root, parent, ppage, 1)
    _ekw_pgwi(h, parent, ppage)
    _ekw_pgwi(h, left, c1page)
    _ekw_pgwi(h, right, c2page)
    _ekw_pgwi(h, new, c3page)
    return overfl
end

"""
    _ekw_trrk(h, tree, left, right, parent, pkidx, nrot) -> Nothing

Rotate `nrot` keys from the node `left` through the parent key number `pkidx` into its
right sibling `right` (`zzektrrk`); a negative `nrot` rotates to the left. Precondition:
both nodes keep 40 to 63 keys, else `SpiceError(:bug)`. The words of a node after its
last key are not cleared, as in CSPICE.
"""
function _ekw_trrk(h::DASHandle, tree::Integer, left::Int, right::Int, parent::Int,
        pkidx::Int, nrot::Int)
    nrot == 0 && return nothing
    root = Int(tree)
    (left == root || right == root) &&
        _ekw_bug("Input node is root; only children are eligible for key rotation.")
    lpage = _ekw_pgri(h, left)
    rpage = _ekw_pgri(h, right)
    ppage = _ekw_pgri(h, parent)
    keybas, datbas, kidbas = _ekw_bases(parent == root)
    _ekw_sibling_check(ppage, kidbas, pkidx, left, right, parent)
    lnkeys = lpage[0]
    rnkeys = rpage[0]
    if nrot > 0
        (lnkeys - nrot < 40 || rnkeys + nrot > 63) && _ekw_bug("Node $left and right sibling " *
            "$right contain $lnkeys and $rnkeys keys respectively; rotation of $nrot keys " *
            "to the right will violate the key count bounds of 40:63.")
    else
        (lnkeys - nrot > 63 || rnkeys + nrot < 40) && _ekw_bug("Node $left and right sibling " *
            "$right contain $lnkeys and $rnkeys keys respectively; rotation of $(-nrot) " *
            "keys to the left will violate the key count bounds of 40:63.")
    end
    lnsize = pkidx == 1 ? ppage[keybas] - 1 : ppage[keybas+pkidx-1] - ppage[keybas+pkidx-2] - 1
    if nrot > 0
        schlep = nrot
        remain = lnkeys - schlep
        futrpk = lpage[remain+1]
        subsiz = lnsize - futrpk
        drotat = -futrpk
        dshift = subsiz + 1
        nvopar = dshift
        dpar = -dshift
        for i in rnkeys:-1:1
            rpage[i+1+schlep-1] = rpage[i] + dshift
        end
        for i in rnkeys:-1:1
            rpage[i+128+schlep-1] = rpage[i+127]
        end
        for i in rnkeys+1:-1:1
            rpage[i+64+schlep-1] = rpage[i+63]
        end
        rpage[schlep] = nvopar
        rpage[schlep+127] = ppage[datbas+pkidx-1]
        ppage[datbas+pkidx-1] = lpage[remain+128]
        ppage[keybas+pkidx-1] = ppage[keybas+pkidx-1] + dpar
        for i in 1:schlep-1
            rpage[i] = lpage[remain+2+i-1] + drotat
        end
        _ekw_move!(rpage, 128, lpage, remain + 129, schlep - 1)
        _ekw_move!(rpage, 64, lpage, remain + 65, schlep)
        lpage[0] = lpage[0] - schlep
        rpage[0] = rpage[0] + schlep
    else
        schlep = -nrot
        remain = rnkeys - schlep
        futrpk = rpage[schlep]
        subsiz = futrpk - 1
        drotat = lnsize + 1
        dshift = -(subsiz + 1)
        nvopar = drotat
        dpar = -dshift
        lpage[lnkeys+1] = nvopar
        lpage[lnkeys+128] = ppage[datbas+pkidx-1]
        ppage[datbas+pkidx-1] = rpage[schlep+127]
        ppage[keybas+pkidx-1] = ppage[keybas+pkidx-1] + dpar
        _ekw_move!(lpage, lnkeys + 2, rpage, 1, schlep - 1)
        _ekw_move!(lpage, lnkeys + 129, rpage, 128, schlep - 1)
        _ekw_move!(lpage, lnkeys + 65, rpage, 64, schlep)
        for i in 1:schlep-1
            lpage[lnkeys+2+i-1] = lpage[lnkeys+2+i-1] + drotat
        end
        for i in 1:remain
            rpage[i] = rpage[i+1+schlep-1] + dshift
        end
        for i in 1:remain
            rpage[i+127] = rpage[i+128+schlep-1]
        end
        for i in 1:remain+1
            rpage[i+63] = rpage[i+64+schlep-1]
        end
        lpage[0] = lpage[0] + schlep
        rpage[0] = rpage[0] - schlep
    end
    _ekw_pgwi(h, parent, ppage)
    _ekw_pgwi(h, left, lpage)
    _ekw_pgwi(h, right, rpage)
    return nothing
end

"""
    _ekw_trbn(h, tree, left, right, parent, pkidx) -> Nothing

Balance the two siblings `left` and `right` so that their key counts differ by at most
one (`zzektrbn`). Precondition: they hold 82 to 124 keys together, else `:bug`.
"""
function _ekw_trbn(h::DASHandle, tree::Integer, left::Int, right::Int, parent::Int, pkidx::Int)
    (left == tree || right == tree) &&
        _ekw_bug("Input node is root; only children can be balanced.")
    lnkeys = _ekw_trnk(h, tree, left)
    rnkeys = _ekw_trnk(h, tree, right)
    sum = lnkeys + rnkeys
    (82 <= sum <= 124) || _ekw_bug("Node $left and right sibling $right contain $lnkeys " *
        "and $rnkeys keys respectively; count sum should be in range 82:124.")
    schlep = lnkeys > rnkeys ? lnkeys - (sum + 1) ÷ 2 :
             lnkeys < rnkeys ? -(rnkeys - (sum + 1) ÷ 2) : 0
    _ekw_trrk(h, tree, left, right, parent, pkidx, schlep)
    return nothing
end

"""
    _ekw_trin(h, tree, key, value) -> Nothing

Insert `value` at the key `key` of the tree with the root page `tree` and restore the
balance (`zzektrin`): an overflow is passed to a sibling that has room, else two
siblings are split into three, else the root is split. Precondition: `1 <= key <=
nkeys + 1` (`:invalidindex`). Postcondition: the tree has one key more; the values of
the old keys `key`, `key + 1`, ... are now at `key + 1`, `key + 2`, ...
"""
function _ekw_trin(h::DASHandle, tree::Integer, key::Integer, value::Integer)
    root = Int(tree)
    lkey = Int(key)
    _ekw_trui(h, tree, lkey, value) || return nothing
    node = left = right = parent = lnode = rnode = pkidx = lpidx = rpidx = 0
    state = 2
    while state != 1
        # invariant: the node that holds `lkey` has one key too many; all other nodes
        # are inside their bounds
        if state == 2
            _, node, _, _, _ = _ek_trlk(h.ints, tree, lkey)
            if node == root
                state = 6
            else
                parent, _, _, lpidx, _, left, rpidx, _, right = _ekw_trpi(h, tree, lkey)
                if left > 0 && _ekw_trnk(h, tree, left) < 62
                    lnode = left
                    rnode = node
                    pkidx = lpidx
                    state = 4
                else
                    state = 3
                end
            end
        elseif state == 3
            if right > 0
                lnode = node
                rnode = right
                pkidx = rpidx
                state = _ekw_trnk(h, tree, right) < 62 ? 4 : 5
            else
                lnode = left
                rnode = node
                pkidx = lpidx
                state = 5
            end
        elseif state == 4
            _ekw_trbn(h, tree, lnode, rnode, parent, pkidx)
            state = 1
        elseif state == 5
            # `trust` is a key that stays in the same child of `parent` after the split.
            trust = node == lnode ? _ekw_trki(h, tree, lkey, 1) :
                    _ekw_trki(h, tree, lkey, _ekw_trnk(h, tree, node))
            if _ekw_tr23(h, tree, lnode, rnode, parent, pkidx)
                if parent == root
                    state = 6
                else
                    parent, pkey, _, lpidx, _, left, rpidx, _, right = _ekw_trpi(h, tree, trust)
                    lkey = pkey
                    state = 2
                end
            else
                state = 1
            end
        else
            _ekw_tr13(h, tree)
            state = 1
        end
    end
    return nothing
end

"Append `value` to the tree `tree`; returns its key (`zzektrap`)."
function _ekw_trap(h::DASHandle, tree::Integer, value::Integer)
    key = _ek_trsz(h.ints, tree) + 1
    _ekw_trin(h, tree, key, value)
    return key
end

# ---------------------------------------------------------------------------------------
# Tree: deletion
# ---------------------------------------------------------------------------------------

"""
    _ekw_trud(h, tree, key) -> (trgkey, underflow)

Delete the key `key` without balancing (`zzektrud`): the keys after it move down by
one. `trgkey` is a key of the node that lost a key (0 if the tree is empty now), and
`underflow` tells that this node holds one key less than its lower bound. A key that
is not in a leaf is replaced by its predecessor, which is taken from a leaf.
Precondition: `1 <= key <= nkeys`, else `:invalidindex`.
"""
function _ekw_trud(h::DASHandle, tree::Integer, key::Int)
    root = Int(tree)
    rpage = _ekw_pgri(h, root)
    totkey = rpage[2]
    (1 <= key <= totkey) || throw(SpiceError(:invalidindex,
        "Key = $key. Valid range is 1:$totkey.  File = $(h.path)."))
    nnode = rpage[1]
    depth = rpage[3]
    if nnode == 1
        nkeys = totkey
        rpage[2] = nkeys - 1
        rpage[4] = nkeys - 1
        for i in key:nkeys-1
            rpage[i+4] = rpage[i+5] - 1
            rpage[i+171] = rpage[i+172]
        end
        for i in key:nkeys
            rpage[i+87] = rpage[i+88]
        end
        rpage[nkeys+4] = 0
        rpage[nkeys+171] = 0
        rpage[nkeys+88] = 0
        nkeys -= 1
        trgkey = nkeys > 0 ? rpage[5] : 0
        _ekw_pgwi(h, root, rpage)
        return trgkey, false
    elseif key == totkey
        _, target, toffst, _, _ = _ek_trlk(h.ints, tree, key)
        tpage = _ekw_pgri(h, target)
        nkeys = tpage[0]
        tpage[nkeys] = 0
        tpage[nkeys+127] = 0
        tpage[nkeys+64] = 0
        tpage[0] = tpage[0] - 1
        rpage[2] = totkey - 1
        trgkey = tpage[1] + toffst
        _ekw_pgwi(h, root, rpage)
        _ekw_pgwi(h, target, tpage)
        return trgkey, nkeys == 41
    end
    _, target, toffst, keyidx, level = _ek_trlk(h.ints, tree, key)
    local undrfl::Bool, trgkey::Int, parent::Int, pkey::Int, rpidx::Int
    tpage = _ekw_Page()
    if level == depth
        parent, pkey, _, _, _, _, rpidx, _, _ = _ekw_trpi(h, tree, key)
        _ekw_pgri!(tpage, h, target)
        tnkeys = tpage[0]
        tnkeys < 41 && _ekw_bug("Node = $target. Tree = $tree. File = $(h.path). Key count " *
            "= $tnkeys; max allowed, including overflow, is 41.")
        for i in keyidx:tnkeys-1
            tpage[i] = tpage[i+1] - 1
        end
        for i in keyidx:tnkeys-1
            tpage[i+127] = tpage[i+128]
        end
        for i in keyidx:tnkeys
            tpage[i+63] = tpage[i+64]
        end
        tpage[0] = tnkeys - 1
        undrfl = tnkeys == 41
        trgkey = tpage[1] + toffst
        _ekw_pgwi(h, target, tpage)
    else
        _, leaf, loffst, prev, _ = _ek_trlk(h.ints, tree, key - 1)
        parent, pkey, _, _, _, _, rpidx, _, _ = _ekw_trpi(h, tree, key - 1)
        lpage = _ekw_pgri(h, leaf)
        if target == root
            rpage[keyidx+171] = lpage[prev+127]
        else
            _ekw_pgri!(tpage, h, target)
            tpage[keyidx+127] = lpage[prev+127]
        end
        nlkeys = lpage[0]
        for i in prev:nlkeys-1
            lpage[i] = lpage[i+1] - 1
            lpage[i+127] = lpage[i+128]
        end
        lpage[0] = nlkeys - 1
        undrfl = nlkeys == 41
        trgkey = lpage[1] + loffst
        _ekw_pgwi(h, leaf, lpage)
        target != root && _ekw_pgwi(h, target, tpage)
    end
    while parent != root
        # invariant: as in `_ekw_trui`, with the keys after the deleted key moving down
        paren2, pkey2, _, _, _, _, rpidx2, _, _ = _ekw_trpi(h, tree, pkey)
        if rpidx > 0
            _ekw_pgri!(tpage, h, parent)
            for i in rpidx:tpage[0]
                tpage[i] = tpage[i] - 1
            end
            _ekw_pgwi(h, parent, tpage)
        end
        parent = paren2
        pkey = pkey2
        rpidx = rpidx2
    end
    if rpidx > 0
        for i in rpidx:rpage[4]
            rpage[i+4] = rpage[i+4] - 1
        end
    end
    rpage[2] = totkey - 1
    _ekw_pgwi(h, root, rpage)
    return trgkey, undrfl
end

"""
    _ekw_tr31(h, tree) -> Nothing

Join the root that holds one key and its two children, which hold 81 keys together,
into a root with 82 keys (`zzektr31`). The children are freed; the depth goes down by
one.
"""
function _ekw_tr31(h::DASHandle, tree::Integer)
    root = Int(tree)
    rpage = _ekw_pgri(h, root)
    nrkeys = rpage[4]
    nrkeys == 1 || _ekw_bug("Number of keys in root = $nrkeys; should be 1.")
    child1 = rpage[88]
    child2 = rpage[89]
    c1page = _ekw_pgri(h, child1)
    c2page = _ekw_pgri(h, child2)
    nlkeys = c1page[0]
    nrkeys = c2page[0]
    nlkeys + nrkeys == 81 || _ekw_bug("Number of keys in nodes LEFT = $nlkeys; in RIGHT = " *
        "$nrkeys; counts summing to 81 were expected.")
    rpage[nlkeys+5] = rpage[5]
    rpage[nlkeys+172] = rpage[172]
    _ekw_move!(rpage, 5, c1page, 1, nlkeys)
    _ekw_move!(rpage, 172, c1page, 128, nlkeys)
    _ekw_move!(rpage, 88, c1page, 64, nlkeys + 1)
    middle = nlkeys + 1
    delta = rpage[middle+4]
    for i in 1:nrkeys
        rpage[middle+5+i-1] = c2page[i] + delta
    end
    _ekw_move!(rpage, middle + 172, c2page, 128, nrkeys)
    _ekw_move!(rpage, middle + 88, c2page, 64, nrkeys + 1)
    rpage[4] = 82
    rpage[3] = rpage[3] - 1
    rpage[1] = rpage[1] - 2
    _ekw_pgwi(h, root, rpage)
    _ekw_pgfr(h, 3, child1)
    _ekw_pgfr(h, 3, child2)
    return nothing
end

"""
    _ekw_tr32(h, tree, left, middle, right, parent, lpkidx) -> underflow

Join the three siblings `left`, `middle` and `right`, which hold 122 keys together,
into `left` (61 keys) and `right` (62 keys), and free `middle` (`zzektr32`). `lpkidx`
is the position in `parent` of the key between `left` and `middle`. The result tells
that the parent now holds one key less than its lower bound (for the root: no key).
"""
function _ekw_tr32(h::DASHandle, tree::Integer, left::Int, middle::Int, right::Int,
        parent::Int, lpkidx::Int)
    c1page = _ekw_pgri(h, left)
    c2page = _ekw_pgri(h, middle)
    c3page = _ekw_pgri(h, right)
    ppage = _ekw_pgri(h, parent)
    root = Int(tree)
    keybas, datbas, kidbas = _ekw_bases(parent == root)
    psize = parent == root ? ppage[4] : ppage[0]
    (1 <= lpkidx <= psize - 1) || _ekw_bug("Left parent key of MIDDLE is out of range.  " *
        "Value is $lpkidx; valid range is 1:$(psize - 1)")
    lpkey = ppage[keybas+lpkidx-1]
    rpkey = ppage[keybas+lpkidx]
    lsib = ppage[kidbas+lpkidx-1]
    msib = ppage[kidbas+lpkidx]
    rsib = ppage[kidbas+lpkidx+1]
    (lsib == left && msib == middle && rsib == right) || _ekw_bug("LEFT, RIGHT, MIDDLE, " *
        "PARENT, and PKIDX are inconsistent. LEFT = $left; MIDDLE = $middle; RIGHT = " *
        "$right; PARENT = $parent; LPKIDX = $lpkidx; LSIB derived from PARENT = $lsib; " *
        "MSIB = $msib; RSIB = $rsib.")
    npkeys = psize
    nlkeys = c1page[0]
    nmkeys = c2page[0]
    nrkeys = c3page[0]
    nlkeys + nmkeys + nrkeys == 122 || _ekw_bug("Number of keys in nodes LEFT = $nlkeys; " *
        "in MIDDLE = $nmkeys; in RIGHT = $nrkeys; counts summing to 122 were expected.")
    leftsz = lpkidx == 1 ? lpkey - 1 : lpkey - ppage[keybas+lpkidx-2] - 1
    n = 61 - (nlkeys + 1) + 1
    lmidsz = c2page[n] - 1
    midsiz = rpkey - lpkey - 1
    rmidsz = midsiz - lmidsz - 1
    c1page[nlkeys+1] = leftsz + 1
    c1page[nlkeys+128] = ppage[datbas+lpkidx-1]
    for i in 1:n-1
        c1page[nlkeys+2+i-1] = c2page[i] + leftsz + 1
    end
    _ekw_move!(c1page, nlkeys + 129, c2page, 128, n - 1)
    _ekw_move!(c1page, nlkeys + 65, c2page, 64, n)
    c1page[0] = 61
    rshift = nmkeys - n + 1
    for i in nrkeys:-1:1
        c3page[rshift+1+i-1] = c3page[i] + rmidsz + 1
    end
    for i in nrkeys:-1:1
        c3page[rshift+128+i-1] = c3page[i+127]
    end
    for i in nrkeys+1:-1:1
        c3page[rshift+64+i-1] = c3page[i+63]
    end
    c3page[rshift] = rmidsz + 1
    c3page[rshift+127] = ppage[datbas+lpkidx]
    for i in 1:rshift-1
        c3page[i] = c2page[n+1+i-1] - lmidsz - 1
    end
    _ekw_move!(c3page, 128, c2page, n + 128, rshift - 1)
    _ekw_move!(c3page, 64, c2page, n + 64, rshift)
    c3page[0] = 62
    ppage[keybas+lpkidx-1] = lpkey + lmidsz + 1
    ppage[datbas+lpkidx-1] = c2page[n+127]
    for i in lpkidx+1:npkeys-1
        ppage[keybas+i-1] = ppage[keybas+i]
    end
    for i in lpkidx+1:npkeys-1
        ppage[datbas+i-1] = ppage[datbas+i]
    end
    for i in lpkidx+1:npkeys
        ppage[kidbas+i-1] = ppage[kidbas+i]
    end
    ppage[keybas+npkeys-1] = 0
    ppage[datbas+npkeys-1] = 0
    ppage[kidbas+npkeys] = 0
    local undrfl::Bool
    if parent == root
        ppage[4] = ppage[4] - 1
        undrfl = ppage[4] == 0
    else
        ppage[0] = ppage[0] - 1
        undrfl = ppage[0] == 40
    end
    _ekw_nnode!(h, root, parent, ppage, -1)
    _ekw_pgwi(h, parent, ppage)
    _ekw_pgwi(h, left, c1page)
    _ekw_pgwi(h, right, c3page)
    _ekw_pgfr(h, 3, middle)
    return undrfl
end

"""
    _ekw_trdl(h, tree, key) -> Nothing

Delete the key `key` of the tree with the root page `tree` and restore the balance
(`zzektrdl`): an underflow takes keys from a sibling that has spare keys (also through
the sibling beyond, for a node at an end), else three siblings are joined into two,
else the two children of the root are joined into the root. Precondition: `1 <= key <=
nkeys` (`:invalidindex`). Postcondition: the tree has one key less; the values of the
old keys `key + 1`, `key + 2`, ... are now at `key`, `key + 1`, ...
"""
function _ekw_trdl(h::DASHandle, tree::Integer, key::Integer)
    trgkey, undrfl = _ekw_trud(h, tree, Int(key))
    undrfl || return nothing
    node = lsib = rsib = lkey = rkey = lnkey = lnode = mnode = rnode = 0
    state = 2
    while state != 1
        # invariant: the node that holds `trgkey` has one key too few; all other nodes
        # are inside their bounds
        if state == 2
            _, node, _, _, _ = _ek_trlk(h.ints, tree, trgkey)
            lsib, lkey, rsib, rkey = _ekw_trsb(h, tree, trgkey)
            if lsib > 0
                if _ekw_trnk(h, tree, lsib) > 41
                    lnkey = lkey
                    lnode = lsib
                    rnode = node
                    state = 4
                elseif rsib > 0
                    state = 3
                else
                    state = 7
                end
            else
                state = 3
            end
        elseif state == 3
            if _ekw_trnk(h, tree, rsib) > 41
                lnkey = trgkey
                lnode = node
                rnode = rsib
                state = 4
            elseif lsib > 0
                lnkey = lkey
                lnode = lsib
                mnode = node
                rnode = rsib
                state = 5
            else
                state = 8
            end
        elseif state == 7
            # The node is the last child: look at the sibling beyond the left sibling.
            llsib, llkey, _, _ = _ekw_trsb(h, tree, lkey)
            if llsib > 0
                if _ekw_trnk(h, tree, llsib) > 41
                    parent, _, _, lpidx, _, llsib, _, _, _ = _ekw_trpi(h, tree, lkey)
                    _ekw_trrk(h, tree, llsib, lsib, parent, lpidx, 1)
                    lnkey = lkey
                    lnode = lsib
                    rnode = node
                    state = 4
                else
                    lnkey = llkey
                    lnode = llsib
                    mnode = lsib
                    rnode = node
                    state = 5
                end
            else
                lnode = lsib
                rnode = node
                state = 6
            end
        elseif state == 8
            # The node is the first child: look at the sibling beyond the right sibling.
            _, _, rrsib, _ = _ekw_trsb(h, tree, rkey)
            if rrsib > 0
                if _ekw_trnk(h, tree, rrsib) > 41
                    parent, _, _, _, _, _, rpidx, _, rrsib = _ekw_trpi(h, tree, rkey)
                    _ekw_trrk(h, tree, rsib, rrsib, parent, rpidx, -1)
                    lnkey = trgkey
                    lnode = node
                    rnode = rsib
                    state = 4
                else
                    lnkey = trgkey
                    lnode = node
                    mnode = rsib
                    rnode = rrsib
                    state = 5
                end
            else
                lnode = node
                rnode = rsib
                state = 6
            end
        elseif state == 4
            parent, _, _, _, _, _, rpidx, _, _ = _ekw_trpi(h, tree, lnkey)
            _ekw_trbn(h, tree, lnode, rnode, parent, rpidx)
            state = 1
        elseif state == 5
            trust = _ekw_trki(h, tree, lnkey, 1)
            parent, _, _, _, _, _, rpidx, _, _ = _ekw_trpi(h, tree, lnkey)
            if _ekw_tr32(h, tree, lnode, mnode, rnode, parent, rpidx)
                _, pkey, _, _, _, _, _, _, _ = _ekw_trpi(h, tree, trust)
                trgkey = pkey
                state = 2
            else
                state = 1
            end
        else
            _ekw_tr31(h, tree)
            state = 1
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Tree: one-shot load and release
# ---------------------------------------------------------------------------------------

"""
    _ekw_tr1s(h, tree, values) -> Nothing

Load the empty tree `tree` with `values` in one step (`zzektr1s`): the keys 1 to
`length(values)` in a tree of the smallest depth, whose subtrees of one node differ in
size by at most one. Precondition: the tree is empty (`:nonemptytree`) and the depth
is at most 10 (`:counttoolarge`).

CSPICE fills every node in one page buffer, and it does not clear the buffer before it
fills a leaf, so a leaf keeps words of the node that the buffer held before (its
parent). The lookups never read those words. The same buffer is used here, so the
files are equal.
"""
function _ekw_tr1s(h::DASHandle, tree::Integer, values::AbstractVector{<:Integer})
    size = length(values)
    tsize = _ek_trsz(h.ints, tree)
    tsize > 0 && throw(SpiceError(:nonemptytree,
        "Tree has size $tsize; should be empty.EK = $(h.path); TREE = $tree."))
    tsize = 82
    d = 1
    s = 1
    while tsize < size                   # invariant: a tree of depth d holds up to tsize keys
        d += 1
        s = d == 2 ? 82 : s * 63 + 83
        tsize = s * 62 + 82
    end
    d > _EK_TRMXDP && throw(SpiceError(:counttoolarge, "Tree has depth $d; max supported " *
        "depth is $_EK_TRMXDP.EK = $(h.path); TREE = $tree."))
    stnode = zeros(Int, _EK_TRMXDP)
    stnkey = zeros(Int, _EK_TRMXDP)
    stlsiz = zeros(Int, _EK_TRMXDP)
    stnbig = zeros(Int, _EK_TRMXDP)
    stnext = zeros(Int, _EK_TRMXDP)
    stnbas = zeros(Int, _EK_TRMXDP)
    page = _ekw_Page()
    node = Int(tree)
    subsiz = size
    next = 1
    level = 1
    basidx = 0
    totnod = 0
    nkeys = nkids = bigsiz = nbig = 0
    while level > 0
        # invariant: the nodes left of the path to `node` are complete; `basidx` keys
        # precede the subtree of `node`, which will hold `subsiz` keys
        if next == 1
            if level == 1
                if d > 1
                    nnodes = 0
                    for _ in 1:d-1
                        nnodes = nnodes * 63 + 1
                    end
                    maxsiz = nnodes * 62
                    q = subsiz + 1
                    div = maxsiz + 1
                    nkids = (q + div - 1) ÷ div
                    n = subsiz - (nkids - 1)
                    bigsiz = (n + nkids - 1) ÷ nkids
                    nsmall = nkids * bigsiz - n
                    nbig = nkids - nsmall
                    nkeys = nkids - 1
                else
                    nkeys = size
                    nkids = 0
                end
                _ekw_pgri!(page, h, tree)
                for i in 1:nkeys
                    npred = i == 1 ? 0 : page[i+3]
                    key = d > 1 ? (i <= nbig ? npred + bigsiz + 1 : npred + bigsiz) : i
                    page[i+4] = key
                    page[i+171] = values[key]
                end
                totnod = 1
                for i in 1:nkids
                    child, _ = _ekw_pgal(h, 3)
                    page[i+87] = child
                    totnod += 1
                end
                page[2] = size
                page[3] = d
                page[4] = nkeys
                page[1] = 0
                _ekw_pgwi(h, tree, page)
            elseif level < d
                _ekw_clear!(page, 0, 256)
                n = subsiz - 41
                reqsiz = (n + 40) ÷ 41
                nnodes = 0
                for _ in 1:d-level
                    nnodes = nnodes * 63 + 1
                end
                maxsiz = nnodes * 62
                if reqsiz <= maxsiz
                    nkids = 42
                    bigsiz = reqsiz
                    nsmall = bigsiz * nkids - n
                    nbig = nkids - nsmall
                else
                    q = subsiz + 1
                    div = maxsiz + 1
                    nkids = (q + div - 1) ÷ div
                    q = subsiz - nkids + 1
                    div = nkids
                    bigsiz = (q + div - 1) ÷ div
                    nsmall = bigsiz * nkids - q
                    nbig = nkids - nsmall
                end
                nkeys = nkids - 1
                for i in 1:nkeys
                    npred = i == 1 ? basidx : basidx + page[i-1]
                    key = i <= nbig ? npred + bigsiz + 1 : npred + bigsiz
                    page[i] = key - basidx
                    page[i+127] = values[key]
                end
                for i in 1:nkids
                    child, _ = _ekw_pgal(h, 3)
                    page[i+63] = child
                    totnod += 1
                end
                page[0] = nkeys
                _ekw_pgwi(h, node, page)
            end
            if level < d
                stnode[level] = node
                stnkey[level] = nkeys
                stlsiz[level] = bigsiz
                stnbig[level] = nbig
                stnext[level] = 2
                stnbas[level] = basidx
                kidbas = level == 1 ? 88 : 64
                level += 1
                node = page[kidbas]
                subsiz = bigsiz
            elseif level > 1
                # A leaf: the buffer still holds the words of the parent.
                nkeys = subsiz
                for i in 1:nkeys
                    page[i] = i
                    page[i+127] = values[basidx+i]
                end
                page[0] = nkeys
                _ekw_pgwi(h, node, page)
                level -= 1
                node = stnode[level]
                nkeys = stnkey[level]
                bigsiz = stlsiz[level]
                nbig = stnbig[level]
                next = stnext[level]
                basidx = stnbas[level]
                nkids = nkeys + 1
                _ekw_pgri!(page, h, node)
            else
                level = 0
            end
        elseif next <= nkids
            stnext[level] = next + 1
            kidbas = level == 1 ? 88 : 64
            node = page[kidbas+next-1]
            subsiz = next <= nbig ? stlsiz[level] : stlsiz[level] - 1
            basidx = if next <= nbig + 1
                stnbas[level] + (next - 1) * stlsiz[level] + (next - 1)
            else
                stnbas[level] + nbig * stlsiz[level] + (next - nbig - 1) * (stlsiz[level] - 1) +
                (next - 1)
            end
            level += 1
            next = 1
        else
            level -= 1
            if level >= 1
                node = stnode[level]
                nkeys = stnkey[level]
                bigsiz = stlsiz[level]
                nbig = stnbig[level]
                next = stnext[level]
                basidx = stnbas[level]
                nkids = nkeys + 1
                _ekw_pgri!(page, h, node)
            end
        end
    end
    h.ints[_ek_ibase(tree)+2] = totnod
    return nothing
end

"""
    _ekw_trfr(h, tree) -> Nothing

Free all pages of the tree with the root page `tree` (`zzektrfr`). Precondition: the
depth is at most 10, else `SpiceError(:invalidformat)`.
"""
function _ekw_trfr(h::DASHandle, tree::Integer)
    page = _ekw_pgri(h, tree)
    depth = page[3]
    depth > _EK_TRMXDP && throw(SpiceError(:invalidformat, "Tree has depth $depth; max " *
        "supported depth is $_EK_TRMXDP.EK = $(h.path); TREE = $tree."))
    stack = zeros(Int, 3 * _EK_TRMXDP)
    level = 1
    remain = page[1]
    node = Int(tree)
    nkids = depth == 1 ? 0 : page[4] + 1
    first = 1
    while remain > 0                     # invariant: `remain` nodes of the tree are not free
        if nkids > 0
            kidbas = level == 1 ? 88 : 64
            stack[3level-2] = node
            stack[3level-1] = nkids
            stack[3level] = first
            level += 1
            node = page[kidbas+first-1]
            _ekw_pgri!(page, h, node)
            nkids = level < depth ? page[0] + 1 : 0
            first = 1
        else
            _ekw_pgfr(h, 3, node)
            remain -= 1
            level -= 1
            if level > 0
                node = stack[3level-2]
                first = stack[3level] + 1
                nkids = stack[3level-1] - 1
                _ekw_pgri!(page, h, node)
            end
        end
    end
    return nothing
end
