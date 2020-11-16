using Gamess
using Test

load_file(name) = Gamess.load(joinpath(dirname(@__FILE__), "data", "$(name).out"))

@testset "Gamess.jl" begin
    @testset "Water CIS" begin
        water = load_file("water")
        atoms = water.molecule.atoms
        @test size(atoms) == (3, 5)
        @test atoms.Atom == ["OXYGEN", "HYDROGEN", "HYDROGEN"]

        ci = water.ci
        states = ci.excited_states
        @test size(states) == (100,3)

        dipoles = ci.length_gauge_dipole
        @test size(dipoles) == (2601,10)
        @test dipoles[1, :z] ≈ 0.958273779741262
        @test dipoles[2599,:z] ≈ 0.044333
        @test dipoles[2599,:A] ≈ -5.4312e-41
        @test dipoles[2599,:B] ≈ 797650.0
    end

    @testset "Neon CIS" begin
        neon = load_file("ne+-guga")
        atoms = neon.molecule.atoms
        @test size(atoms) == (1, 5)
        @test atoms.Atom == ["NE"]

        ci = neon.ci
        states = ci.excited_states
        @test size(states) == (5,2)

        dipoles = ci.length_gauge_dipole
        @test size(dipoles) == (10,10)
        @test dipoles[8, :z] ≈ -0.521326

        dipoles = ci.velocity_gauge_dipole
        @test size(dipoles) == (10,8)
    end
end
