
struct Component
    variable::Int64
    poweridx::Int64
    s::Float64
    c::Float64
end

struct Frequency
    coefficients_idx::Array{Int64, 1}
    coefficients::SVector{14, Float64}
    components::Array{Component, 1}
end

function assemble_frequencies(tables)

    frequencies = Array{Frequency, 1}(0)

    for k in 1:length(tables)

        for i in 1:size(tables[k], 1)
            row = tables[k][i, :]

            comp = Component(k, row[1]+1, row[3], row[4])

            j = findfirst(f -> f.coefficients == row[5:18], frequencies)

            if j !== 0
                push!(frequencies[j].components, comp)
            else
                freq = Frequency(find(!iszero, row[5:18]), row[5:18], [comp])
                push!(frequencies, freq)
            end
        end

    end

    return SVector{length(frequencies), Frequency}(frequencies)
end

data_X = sortrows(readdlm(joinpath(@__DIR__, "data/nutation_2000A_X.dlm")), by = r -> hypot(r[3], r[4]))
data_Y = sortrows(readdlm(joinpath(@__DIR__, "data/nutation_2000A_Y.dlm")), by = r -> hypot(r[3], r[4]))
data_s = sortrows(readdlm(joinpath(@__DIR__, "data/nutation_2000A_s.dlm")), by = r -> hypot(r[3], r[4]))

const frequencies = assemble_frequencies([data_X, data_Y, data_s])
