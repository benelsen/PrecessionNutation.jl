using Base.Test
using BenchmarkTools
using AstronomicalTime

# using Revise

using PrecessionNutation

const TENTH_MICROARCSEC = deg2rad(0.1 / 3600e6)

@printf "Testing accuracy magnitude: %6.2f (0.1 μas)\n" log10(TENTH_MICROARCSEC)

@testset "PrecessionNutation" begin

    @testset "Validating PrecessionNutation against ERFA" begin

        @testset "IAU 2006 precession and IAU 2000A_R06 nutation CIO XYs" begin

            @testset "$year" for year in [2020] # [1980, 2000, 2020, 2050, 2100]
                ep = TTEpoch("$(year)-01-01T12:00:00.000")

                x, y, s = [precession_nutation_06(ep)...]
                x_erfa, y_erfa, s_erfa = [precession_nutation_erfa(ep)...]

                @printf "Differences in Year %4d: \n" year
                @printf "       PN.jl -      ERFA =      Diff | Rel. Err. | <= 0.1 μas \n"
                @printf "x:  1.91e-03 -  1.91e-03 =  2.17e-19 |     0.000 | true \n"
                @printf "X: %9.2e - %9.2e = %9.2e | %9.3f | %s \n" x x_erfa (x-x_erfa) (x-x_erfa)/x_erfa norm(x-x_erfa) <= TENTH_MICROARCSEC
                @printf "Y: %9.2e - %9.2e = %9.2e | %9.3f | %s \n" y y_erfa (y-y_erfa) (y-y_erfa)/y_erfa norm(y-y_erfa) <= TENTH_MICROARCSEC
                @printf "s: %9.2e - %9.2e = %9.2e | %9.3f | %s \n" s s_erfa (s-s_erfa) (s-s_erfa)/s_erfa norm(s-s_erfa) <= TENTH_MICROARCSEC

                @test x ≈ x_erfa rtol=0 atol=TENTH_MICROARCSEC
                @test y ≈ y_erfa rtol=0 atol=TENTH_MICROARCSEC
                @test s ≈ s_erfa rtol=0 atol=TENTH_MICROARCSEC
            end

        end

    end

    @testset "Benchmarking against ERFA" begin
        ep = TTEpoch("2018-04-17T18:31:20.921")

        println("\n", repeat("=", 20))
        println("ERFA: \n")
        markERFA = @benchmark precession_nutation_erfa($ep)

        display(markERFA)
        println("\n")

        println("\n", repeat("=", 20))
        println("Variant 1: \n")

        markPN = @benchmark precession_nutation_06($ep)
        display(markPN)
        println("\n")

        judgement = judge(median(markPN), median(markERFA), time_tolerance = 0.5, memory_tolerance = 5.00)
        display(judgement)
        println("\n")

        @test !isregression(judgement)

    end
end
