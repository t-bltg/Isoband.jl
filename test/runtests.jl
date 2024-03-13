using Isoband
using Test

# let's use this quadruple bar thing as the replacement for setequal() from the R tests for convenience
≣(s1::Set, s2::Set) = s1 == s2
≣(s1, s2) = Set(s1) == Set(s2)

"float range"
fr(r::UnitRange) = range(start=F(r.start), stop=F(r.stop))
fr(r::StepRange) = range(start=F(r.start), stop=F(r.stop), step=F(r.step))

F = Float32
include("float_tests.jl")

F = Float64
include("float_tests.jl")
