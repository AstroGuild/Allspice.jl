# commnt: read, add and delete the comments of a binary kernel.

"Architecture of a binary kernel, `:daf` or `:das`; throws for any other file."
function _binarch(path::AbstractString)
    arch, _ = getfat(path)
    arch == "DAF" && return :daf
    arch == "DAS" && return :das
    throw(ArgumentError("$path is not a binary kernel (DAF or DAS file); it is $arch"))
end

"""
    comments(path) -> Vector{String}

Return the lines of the comment area of the binary kernel `path`, as `commnt -r`. The
file can be a DAF file (SPK, CK, binary PCK) or a DAS file (DSK, EK). The result is
empty if the file has no comments.

Throws `ArgumentError` if `path` is not a DAF or DAS file.
"""
function comments(path::AbstractString)
    if _binarch(path) === :daf
        h = dafopr(path)
        try
            return dafec(h)
        finally
            dafcls(h)
        end
    else
        h = dasopr(path)
        try
            return dasec(h)
        finally
            dascls(h)
        end
    end
end

"""
    addcomments(path, lines) -> Nothing

Append `lines` to the comment area of the binary kernel `path`, as `commnt -a`. The
file is changed in place. Each line holds printable ASCII characters only.

Throws `ArgumentError` if `path` is not a DAF or DAS file, and the `SpiceError` of
`Allspice.dafac` or `Allspice.dasac` for an empty `lines` or a character that is not
printable.
"""
function addcomments(path::AbstractString, lines::AbstractVector{<:AbstractString})
    if _binarch(path) === :daf
        h = dafopw(path)
        try
            dafac(h, lines)
        finally
            dafcls(h)
        end
    else
        h = dasopw(path)
        try
            dasac(h, lines)
        finally
            dascls(h)
        end
    end
    return nothing
end

"""
    deletecomments(path) -> Nothing

Delete the comment area of the binary kernel `path`, as `commnt -d`. The file is changed
in place.

Throws `ArgumentError` if `path` is not a DAF or DAS file.
"""
function deletecomments(path::AbstractString)
    if _binarch(path) === :daf
        h = dafopw(path)
        try
            dafdc(h)
        finally
            dafcls(h)
        end
    else
        h = dasopw(path)
        try
            dasdc(h)
        finally
            dascls(h)
        end
    end
    @assert isempty(comments(path))
    return nothing
end
