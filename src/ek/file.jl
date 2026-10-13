# The EK file layer: pages, trees, segment and column descriptors.
#
# An EK file of the current architecture (code 8, "paged DAS") divides each of the three
# DAS address spaces into pages: 1024 characters, 128 doubles, 256 integers. The first
# 256 integers hold the paging metadata, so integer page `p` starts after the address
# `256 p`. The last words of a data page hold a forward pointer and a link count:
#
#   character page   characters 1 to 1014 data, 1015 to 1019 forward pointer, 1020 to 1024
#                    link count, both as five base-128 digits
#   double page      doubles 1 to 126 data, 127 forward pointer, 128 link count
#   integer page     integers 1 to 254 data, 255 forward pointer, 256 link count
#
# Integer page 1 is the file metadata: its first integer is the root page of the segment
# tree. A tree is a B*-tree whose nodes are integer pages; it maps the keys 1 to n to
# integers. The segment tree maps a segment number to the integer page that starts with
# the segment descriptor (24 integers) and the column descriptors (11 integers each).
#
# The functions here read from a "store": any object with the fields `chars`, `dps` and
# `ints` (a `DASHandle` or an `EKFile`). They hold no state, so many tasks can read one
# store at the same time.

const _EK_ARCH = 8                # architecture code at integer address 1
const _EK_SDSCSZ = 24             # size of a segment descriptor
const _EK_CDSCSZ = 11             # size of a column descriptor
const _EK_TNAMSZ = 64             # characters of a table name
const _EK_CNAMSZ = 32             # characters of a column name
const _EK_MXCLSG = 100            # largest number of columns of a segment
const _EK_CPSIZE = 1014           # data characters of a character page
const _EK_DPSIZE = 126            # data doubles of a double page
const _EK_IPSIZE = 254            # data integers of an integer page
const _EK_NULL = -2               # data pointer of a null entry
const _EK_UNINIT = -1             # data pointer of an entry that was never written
const _EK_NOBACK = -3             # data pointer without backup (shadowing; never written here)

# Positions in a segment descriptor (`eksegdsc.inc`).
const _EK_EKTIDX = 1              # segment type
const _EK_SNOIDX = 2              # segment number
const _EK_IMDIDX = 3              # integer base of the segment metadata
const _EK_TNMIDX = 4              # character base of the table name
const _EK_NCIDX = 5               # number of columns
const _EK_NRIDX = 6               # number of rows
const _EK_RTIDX = 7               # root page of the record tree
const _EK_CPTIDX = 8              # root page of the character data page tree
const _EK_DPTIDX = 9              # root page of the double data page tree
const _EK_IPTIDX = 10             # root page of the integer data page tree
const _EK_MFLIDX = 11             # "modified" flag
const _EK_IFLIDX = 12             # "initialized" flag
const _EK_SHDIDX = 13             # shadowing flag
const _EK_CFHIDX = 14             # companion file handle
const _EK_CSNIDX = 15             # companion segment number
const _EK_LCPIDX = 16             # last character, double, integer data page: 16, 17, 18
const _EK_LCWIDX = 19             # last word in use of those pages: 19, 20, 21
const _EK_NMLIDX = 22             # character base of the column names

# Positions in a column descriptor (`ekcoldsc.inc`).
const _EK_CLSIDX = 1
const _EK_TYPIDX = 2
const _EK_LENIDX = 3
const _EK_SIZIDX = 4
const _EK_NAMIDX = 5
const _EK_IXTIDX = 6
const _EK_IXPIDX = 7
const _EK_NFLIDX = 8
const _EK_ORDIDX = 9
const _EK_METIDX = 10

# Positions in a tree page (`ektree.inc`). The root node holds up to 82 keys, another
# node up to 62; a key of a node is the number of keys of the subtree of the node that
# are not greater than it, so a key is the offset of the node plus the stored value.
const _EK_TRTYPE = 1              # root: tree type (1)
const _EK_TRNNOD = 2              # root: number of nodes
const _EK_TRNKEY = 3              # root: number of keys of the tree
const _EK_TRDPTH = 4              # root: depth
const _EK_TRNKR = 5               # root: number of keys of the root node
const _EK_TRKEYR = 5              # root: base of the keys
const _EK_TRKIDR = 88             # root: base of the child pointers
const _EK_TRDATR = 172            # root: base of the data values
const _EK_TRNKC = 1               # child: number of keys
const _EK_TRKEYC = 1              # child: base of the keys
const _EK_TRKIDC = 64             # child: base of the child pointers
const _EK_TRDATC = 128            # child: base of the data values
const _EK_MXKEYR = 82
const _EK_MXKIDR = 83
const _EK_MXKEYC = 62
const _EK_MXKIDC = 63
const _EK_MNKEYC = 41             # smallest number of keys of a node that is not the root
const _EK_MNKIDC = 42
const _EK_MNKEYR = 1
const _EK_MNKIDR = 2
const _EK_TRMXDP = 10             # largest depth of a tree

@inline _ek_cbase(p::Integer) = (Int(p) - 1) * 1024
@inline _ek_dbase(p::Integer) = (Int(p) - 1) * 128
@inline _ek_ibase(p::Integer) = (Int(p) - 1) * 256 + 256

"The integer that five base-128 characters at the character address `addr` encode (`prtdec`)."
@inline function _ek_gei(chars::AbstractVector{UInt8}, addr::Integer)
    n = 0
    for k in 0:4                         # invariant: n holds the digits before k
        n = 128 * n + Int(chars[addr+k])
    end
    return n
end

"""
    _ek_check(s, path, who) -> Nothing

The check of `zzekpgch`: the store `s` holds a paged EK of architecture 8 whose last
addresses lie inside its allocated pages. Throws `SpiceError(:wrongarchitecture)` or
`SpiceError(:invalidformat)`.
"""
function _ek_check(s, path::AbstractString, who::AbstractString)
    ints = s.ints
    (length(ints) >= 16 && ints[1] == _EK_ARCH) || throw(SpiceError(:wrongarchitecture,
        "$who: file $path has architecture $(isempty(ints) ? 0 : ints[1]), which is invalid " *
        "for paged access; the EK software version is $_EK_ARCH"))
    topc = Int(ints[4]) * 1024
    topd = Int(ints[9]) * 128
    topi = Int(ints[14]) * 256 + 256
    length(s.chars) > topc && throw(SpiceError(:invalidformat,
        "$who: file $path has last char address $(length(s.chars)); `top' = $topc."))
    length(s.dps) > topd && throw(SpiceError(:invalidformat,
        "$who: file $path has last d.p. address $(length(s.dps)); `top' = $topd."))
    length(ints) > topi && throw(SpiceError(:invalidformat,
        "$who: file $path has last int. address $(length(ints)); `top' = $topi."))
    return nothing
end

"Number of elements of the sorted integers `ints[first:first+n-1]` that are not greater than `key` (`lstlei`)."
@inline function _ek_lstlei(ints::AbstractVector{Int32}, first::Int, n::Int, key::Int)
    lo = 0
    hi = n
    while lo < hi                        # invariant: elements 1..lo are <= key, elements hi+1..n are > key
        mid = (lo + hi + 1) >> 1
        if ints[first+mid-1] <= key
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

"Number of keys of the tree with the root page `tree` (`zzektrsz`)."
@inline _ek_trsz(ints::AbstractVector{Int32}, tree::Integer) = Int(ints[_ek_ibase(tree)+_EK_TRNKEY])

"""
    _ek_trlk(ints, tree, key) -> (value, node, noffst, idx, level)

Look up `key` in the tree with the root page `tree`, as `zzektrlk`: the data value, the
page of the node that holds the key, the key offset of that node, the position of the
key in the node, and the level of the node (1 for the root). Precondition:
`1 <= key <= _ek_trsz(ints, tree)`, else `SpiceError(:indexoutofrange)`. Throws
`SpiceError(:bug)` for a damaged tree.
"""
function _ek_trlk(ints::AbstractVector{Int32}, tree::Integer, key::Integer)
    base = _ek_ibase(tree)
    totkey = Int(ints[base+_EK_TRNKEY])
    depth = Int(ints[base+_EK_TRDPTH])
    k = Int(key)
    (1 <= k <= totkey) || throw(SpiceError(:indexoutofrange,
        "Key = $k; valid range = 1:$totkey. Tree = $tree"))
    prev = _ek_lstlei(ints, base + _EK_TRKEYR + 1, Int(ints[base+_EK_TRNKR]), k)
    prvkey = prev > 0 ? Int(ints[base+_EK_TRKEYR+prev]) : 0
    prvkey == k && return Int(ints[base+_EK_TRDATR+prev]), Int(tree), 0, prev, 1
    child = Int(ints[base+_EK_TRKIDR+prev+1])
    noffst = prvkey
    level = 1
    while child > 0
        # invariant: the key is in the subtree of `child`, whose keys start after `noffst`
        level += 1
        level > depth && throw(SpiceError(:bug,
            "Runaway node pointer chain.  Key = $k; valid range = 1:$totkey. Tree = $tree"))
        base = _ek_ibase(child)
        newkey = k - noffst
        prev = _ek_lstlei(ints, base + _EK_TRKEYC + 1, Int(ints[base+_EK_TRNKC]), newkey)
        prvkey = prev > 0 ? Int(ints[base+_EK_TRKEYC+prev]) : 0
        prvkey == newkey && return Int(ints[base+_EK_TRDATC+prev]), child, noffst, prev, level
        child = Int(ints[base+_EK_TRKIDC+prev+1])
        noffst += prvkey
    end
    throw(SpiceError(:bug, "Key $k; valid range = 1:$totkey. Tree = $tree.  Key was not " *
        "found.  This probably indicates a corrupted file or a bug in the EK code."))
end

"Data value of `key` in the tree with the root page `tree` (`zzektrdp`)."
@inline _ek_trdp(ints::AbstractVector{Int32}, tree::Integer, key::Integer) =
    _ek_trlk(ints, tree, key)[1]

"""
    _ek_tree_values(ints, tree) -> Vector{Int32}

All data values of the tree with the root page `tree`, in key order. Postcondition: the
length is the number of keys of the tree. Throws `SpiceError(:bug)` for a tree whose
nodes do not hold the number of keys that the root states.
"""
function _ek_tree_values(ints::AbstractVector{Int32}, tree::Integer)
    n = _ek_trsz(ints, tree)
    out = Vector{Int32}(undef, n)
    depth = Int(ints[_ek_ibase(tree)+_EK_TRDPTH])
    filled = _ek_tree_walk!(out, 0, ints, Int(tree), true, depth)
    filled == n || throw(SpiceError(:bug,
        "the tree with root page $tree holds $filled keys, but its root states $n"))
    return out
end

"""
In-order walk of the node `node`: appends its values to `out` after position `pos` and
returns the new position. `below` is the number of tree levels from this node down (the
depth of the tree for the root); a node of the last level is a leaf, and its child
pointers are not read, because CSPICE can leave old values there.
"""
function _ek_tree_walk!(out::Vector{Int32}, pos::Int, ints::AbstractVector{Int32}, node::Int,
        isroot::Bool, below::Int)
    below > _EK_TRMXDP && throw(SpiceError(:bug, "the tree node chain is deeper than $_EK_TRMXDP"))
    base = _ek_ibase(node)
    nkeys = Int(ints[base+(isroot ? _EK_TRNKR : _EK_TRNKC)])
    kidbase = base + (isroot ? _EK_TRKIDR : _EK_TRKIDC)
    datbase = base + (isroot ? _EK_TRDATR : _EK_TRDATC)
    for i in 1:nkeys                     # invariant: `pos` values are in `out`, in key order
        child = below > 1 ? Int(ints[kidbase+i]) : 0
        child > 0 && (pos = _ek_tree_walk!(out, pos, ints, child, false, below - 1))
        pos < length(out) || throw(SpiceError(:bug, "a tree holds more keys than its root states"))
        pos += 1
        out[pos] = ints[datbase+i]
    end
    child = below > 1 ? Int(ints[kidbase+nkeys+1]) : 0
    child > 0 && (pos = _ek_tree_walk!(out, pos, ints, child, false, below - 1))
    return pos
end

@inline _ek_segtree(s) = Int(s.ints[_ek_ibase(1)+1])

@inline _ek_nseg(s) = _ek_trsz(s.ints, _ek_segtree(s))

"""
    _ek_mloc(s, segno) -> (page, base)

Integer page and integer base address of the metadata of segment `segno` (`zzekmloc`).
Throws `SpiceError(:invalidindex)` for a segment number outside `1:nseg`.
"""
function _ek_mloc(s, segno::Integer)
    nseg = _ek_nseg(s)
    (1 <= segno <= nseg) || throw(SpiceError(:invalidindex,
        "Segment number = $segno; valid range is 1:$nseg."))
    page = _ek_trdp(s.ints, _ek_segtree(s), segno)
    return page, _ek_ibase(page)
end

"The 24 integers of the descriptor of segment `segno`, as a view of `s.ints` (`zzeksdsc`)."
@inline function _ek_sdsc(s, segno::Integer)
    _, base = _ek_mloc(s, segno)
    return view(s.ints, base+1:base+_EK_SDSCSZ)
end

function _ek_rstr(chars::AbstractVector{UInt8}, first::Int, last::Int)
    while last >= first && chars[last] == UInt8(' ')   # invariant: characters after `last` are blank
        last -= 1
    end
    return String(chars[first:last])
end

function _ek_column(s, mbase::Integer, i::Integer)
    b = Int(mbase) + _EK_SDSCSZ + (Int(i) - 1) * _EK_CDSCSZ
    ints = s.ints
    nambas = Int(ints[b+_EK_NAMIDX])
    return _EKColumn(_ek_rstr(s.chars, nambas + 1, nambas + _EK_CNAMSZ), ints[b+_EK_CLSIDX],
        ints[b+_EK_TYPIDX], ints[b+_EK_LENIDX], ints[b+_EK_SIZIDX], ints[b+_EK_IXTIDX],
        ints[b+_EK_IXPIDX], ints[b+_EK_NFLIDX] == 1, ints[b+_EK_ORDIDX], ints[b+_EK_METIDX])
end

"True if `a` and `b` are equal without regard to case and blanks (`eqstr`)."
function _ek_eqstr(a::AbstractString, b::AbstractString)
    ca = codeunits(a)
    cb = codeunits(b)
    i = 1
    j = 1
    up(c) = (UInt8('a') <= c <= UInt8('z')) ? c - 0x20 : c
    while true                           # invariant: the non-blank characters before i and j agree
        while i <= length(ca) && ca[i] == UInt8(' ')
            i += 1
        end
        while j <= length(cb) && cb[j] == UInt8(' ')
            j += 1
        end
        (i > length(ca) || j > length(cb)) && return i > length(ca) && j > length(cb)
        up(ca[i]) == up(cb[j]) || return false
        i += 1
        j += 1
    end
end

"""
    _ek_find_column(s, segdsc, column) -> _EKColumn

The descriptor of the column `column` of the segment with the descriptor `segdsc`, as
`zzekcdsc` (the name is compared without regard to case and blanks). Throws
`SpiceError(:bug)` if the segment has no such column, which is what CSPICE signals.
"""
function _ek_find_column(s, segdsc::AbstractVector{Int32}, column::AbstractString)
    mbase = Int(segdsc[_EK_IMDIDX])
    for i in 1:Int(segdsc[_EK_NCIDX])    # invariant: columns before i do not have the name
        col = _ek_column(s, mbase, i)
        _ek_eqstr(col.name, column) && return col
    end
    throw(SpiceError(:bug, "Descriptor for column $column was not found. " *
        "Segment base = $mbase."))
end

"""
    _ek_segment(s, segno) -> _EKSegment

Decode segment `segno` of the store `s`: descriptor, table name, column descriptors and
the record pointers of all rows. Throws `SpiceError(:invalidindex)` for a bad segment
number, `:invalidtype` for a segment type other than 1 and 2, and `:invalidformat` for
a descriptor whose counts are not possible.
"""
function _ek_segment(s, segno::Integer)
    _, base = _ek_mloc(s, segno)
    segdsc = Vector{Int32}(s.ints[base+1:base+_EK_SDSCSZ])
    ncols = Int(segdsc[_EK_NCIDX])
    nrows = Int(segdsc[_EK_NRIDX])
    stype = Int(segdsc[_EK_EKTIDX])
    (1 <= ncols <= _EK_MXCLSG && nrows >= 0) || throw(SpiceError(:invalidformat,
        "segment $segno has $ncols columns and $nrows rows"))
    tnb = Int(segdsc[_EK_TNMIDX])
    table = _ek_rstr(s.chars, tnb + 1, tnb + _EK_TNAMSZ)
    columns = _EKColumn[_ek_column(s, base, i) for i in 1:ncols]
    recptrs = if stype == 1
        v = _ek_tree_values(s.ints, segdsc[_EK_RTIDX])
        length(v) == nrows || throw(SpiceError(:invalidformat, "segment $segno states " *
            "$nrows rows, but its record tree has $(length(v)) keys"))
        v
    elseif stype == 2
        collect(Int32(1):Int32(nrows))
    else
        throw(SpiceError(:invalidtype, "The segment type $stype is not supported."))
    end
    return _EKSegment(Int32(segno), Int32(stype), table, Int32(nrows), columns, recptrs, segdsc)
end

"Record pointer of row `recno` of the segment with the descriptor `segdsc` (`zzekrplk`)."
function _ek_rplk(s, segdsc::AbstractVector{Int32}, recno::Integer)
    stype = segdsc[_EK_EKTIDX]
    stype == 1 && return _ek_trdp(s.ints, segdsc[_EK_RTIDX], recno)
    stype == 2 || throw(SpiceError(:invalidtype, "The segment type $stype is not supported."))
    return Int(recno)
end

# ---------------------------------------------------------------------------------------
# Loading: the EK part of a `KernelSet`
# ---------------------------------------------------------------------------------------

"""
    open_ek(path) -> EKFile

Open the EK file at `path` for `furnsh` and decode its segments, as the first part of
`eklef`. Contract: `path` names a DAS file that holds a paged EK of architecture 8 with
at least one segment, of either byte order. In the result the address spaces are
memory-mapped where the file allows it. Throws `SpiceError` with the kinds of `dasopr`,
`:wrongarchitecture`, `:invalidformat`, `:eknosegments`, and `:ekcolnummismatch` for a
segment with two columns of the same name.
"""
function open_ek(path::AbstractString)
    h = _das_load(path, "furnsh", false)
    _ek_check(h, path, "furnsh")
    nseg = _ek_nseg(h)
    nseg < 1 && throw(SpiceError(:eknosegments, "The EK file $path contains no segments."))
    segs = _EKSegment[_ek_segment(h, i) for i in 1:nseg]
    for seg in segs
        allunique(c.name for c in seg.columns) || throw(SpiceError(:ekcolnummismatch,
            "The EK file $path could not be loaded; the segment $(seg.segno) contains " *
            "duplicate column names in table $(seg.table)."))
    end
    return EKFile(String(path), :EK, h.chars, h.dps, h.ints, segs)
end
