
@inline as_to_rad(x) = deg2rad(x / 3600.0)
@inline as_to_rad_mod(x) = deg2rad(rem(x, 1296000.0) / 3600.0)
@inline rad_to_as(x) = rad2deg(x) * 3600.0

@inline μas_to_rad(x) = deg2rad(x / 3600e6)
@inline rad_to_μas(x) = rad2deg(x) * 3600e6

@inline jc(ep) = ((julian1(ep) - J2000) + julian2(ep)) / DAYS_PER_CENTURY

function fundamental_arguments(tdb_jc)

    # Mean Anomaly of the Moon (l, F1)
    M_m  = as_to_rad( @evalpoly(tdb_jc,   485_868.249_036, 1_717_915_923.217_8,  31.879_2,  0.051_635, -0.000_244_70) )

    # Mean Anomaly of the Sun (l', F2)
    M_s  = as_to_rad( @evalpoly(tdb_jc, 1_287_104.793_048,   129_596_581.048_1,  -0.553_2,  0.000_136, -0.000_011_49) )

    # Mean Argument of Latitude of the Moon (F, F3)
    u_sm = as_to_rad( @evalpoly(tdb_jc,   335_779.526_232, 1_739_527_262.847_8, -12.751_2, -0.001_037,  0.000_004_17) )

    # Mean Elongation of the Moon from the Sun (D, F4)
    D_s  = as_to_rad( @evalpoly(tdb_jc, 1_072_260.703_692, 1_602_961_601.209_0,  -6.370_6,  0.006_593, -0.000_031_69) )

    # Mean Longitude of the ascending Node of the Moon (Ω, F5)
    Ω_m  = as_to_rad( @evalpoly(tdb_jc,   450_160.398_036,    -6_962_890.543_1,   7.472_2,  0.007_702, -0.000_059_39) )

    # Mean Heliocentric Longitudes of the Planets
    # Souchay et al. (1999)
    λ_M1 = @evalpoly(tdb_jc, 4.402_608_842, 2608.790_314_157_4) # Mercury
    λ_M2 = @evalpoly(tdb_jc, 3.176_146_697, 1021.328_554_621_1) # Venus
    λ_M3 = @evalpoly(tdb_jc, 1.753_470_314,  628.307_584_999_1) # Earth
    λ_M4 = @evalpoly(tdb_jc, 6.203_480_913,  334.061_242_670_0) # Mars
    λ_M5 = @evalpoly(tdb_jc, 0.599_546_497,   52.969_096_264_1) # Jupiter
    λ_M6 = @evalpoly(tdb_jc, 0.874_016_757,   21.329_910_496_0) # Saturn
    λ_M7 = @evalpoly(tdb_jc, 5.481_293_872,    7.478_159_856_7) # Uranus
    λ_M8 = @evalpoly(tdb_jc, 5.311_886_287,    3.813_303_563_8) # Neptune

    # General Precession
    # Kinoshita and Souchay (1990)
    p_λ  = @evalpoly(tdb_jc, 0.0, 0.024_381_75, 0.000_005_386_91)

    return SVector(M_m, M_s, u_sm, D_s, Ω_m, λ_M1, λ_M2, λ_M3, λ_M4, λ_M5, λ_M6, λ_M7, λ_M8, p_λ)
end
