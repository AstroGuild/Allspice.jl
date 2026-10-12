# `@apart` for the test groups: a block that is compiled on its own.
#
# The body of a group is one function, and Julia compiles it in one piece before the
# first test runs. For a body of several hundred lines that takes many seconds. With
# `@apart` before each `@testset` of the body, each test set is a function of its own
# that is compiled when the run reaches it. The measured gain is 10 to 20 percent of the
# run time of the large groups.

"""
    @apart expr

Evaluate `expr` in a function of its own that is called through `Base.invokelatest`, so
that the compiler does not compile `expr` together with the code around it. The value
is the value of `expr`. A local variable of the enclosing scope that `expr` assigns is
shared with it, as in a closure. Do not use it around `return`, `break` or `continue`
that belong to the enclosing code.
"""
macro apart(ex)
    return :(Base.invokelatest(() -> $(esc(ex))))
end
