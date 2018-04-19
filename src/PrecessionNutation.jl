module PrecessionNutation

using AstronomicalTime
using ERFA
using StaticArrays
using MuladdMacro

import EarthOrientation: precession_nutation00

include("helper.jl")
include("tables.jl")

export precession_nutation_erfa, precession_nutation_06

function precession_nutation_erfa(ep::TTEpoch)
    x, y = ERFA.xy06(julian1(ep), julian2(ep))
    s = ERFA.s06(julian1(ep), julian2(ep), x, y)
    x, y, s
end

### IAU 2006 precession and IAU 2000A_R06 nutation
function precession_nutation_IAU2006_IAU2000A_R06(ep_tt::TTEpoch; revision = 0)
    tt_jc = jc(ep_tt)

    if revision === 0
        # P03
        X = @evalpoly(tt_jc, -16_617.0, +2_004_191_898.00,    -429_782.90, -198_618.34,     +7.578,  +5.928_5) # μas
        Y = @evalpoly(tt_jc,  -6_951.0,        -25_896.00, -22_407_274.70,   +1_900.59, +1_112.526,  +0.135_8) # μas
        s = @evalpoly(tt_jc,     +94.0,         +3_808.65,        -122.68,  -72_574.11,    +27.980, +15.620_0) # μas

    elseif revision === 1
        # P03_rev1
        X = @evalpoly(tt_jc, -16_617.0, +2_004_191_804.00,    -429_755.80, -198_618.29,     +7.575,  +5.928_5) # μas
        Y = @evalpoly(tt_jc,  -6_951.0,        -24_867.00, -22_407_272.70,   +1_900.26, +1_112.525,  +0.135_8) # μas

        # TODO: confirm values for s
        s = @evalpoly(tt_jc,     +94.0,         +3_808.65,        -122.68,  -72_574.11,    +27.980, +15.620_0) # μas

    elseif revision === 2
        # P03_rev2
        X = @evalpoly(tt_jc, -16_617.0, +2_004_192_130.00,    -429_775.20, -198_618.39,     +7.576,  +5.928_5) # μas
        Y = @evalpoly(tt_jc,  -6_951.0,        -25_817.00, -22_407_280.10,   +1_900.46, +1_112.526,  +0.135_8) # μas

        # TODO: confirm values for s
        s = @evalpoly(tt_jc,     +94.0,         +3_808.65,        -122.68,  -72_574.11,    +27.980, +15.620_0) # μas

    else
        error("Revision $revision not implemented")
    end

    fund_args = fundamental_arguments(tt_jc)

    X_harm_coeff = zeros(MVector{5, Float64})
    Y_harm_coeff = zeros(MVector{5, Float64})
    s_harm_coeff = zeros(MVector{5, Float64})

    @inbounds for frequency in frequencies

        arg = 0.0
        @inbounds @simd for k in frequency.coefficients_idx
            @muladd arg += frequency.coefficients[k] * fund_args[k]
        end
        si = sin(arg)
        co = cos(arg)

        @inbounds for component in frequency.components
            @muladd ampl = component.s * si + component.c * co

            if component.variable == 1
                X_harm_coeff[component.poweridx] += ampl
            elseif component.variable == 2
                Y_harm_coeff[component.poweridx] += ampl
            else
                s_harm_coeff[component.poweridx] += ampl
            end
        end
    end

    X += @evalpoly(tt_jc, X_harm_coeff[1], X_harm_coeff[2], X_harm_coeff[3], X_harm_coeff[4], X_harm_coeff[5])
    Y += @evalpoly(tt_jc, Y_harm_coeff[1], Y_harm_coeff[2], Y_harm_coeff[3], Y_harm_coeff[4], Y_harm_coeff[5])
    s += @evalpoly(tt_jc, s_harm_coeff[1], s_harm_coeff[2], s_harm_coeff[3], s_harm_coeff[4], s_harm_coeff[5])

    X = μas_to_rad(X)
    Y = μas_to_rad(Y)
    s = μas_to_rad(s)

    X, Y, s - X * Y / 2
end

const precession_nutation_06 = precession_nutation_IAU2006_IAU2000A_R06

end # module
