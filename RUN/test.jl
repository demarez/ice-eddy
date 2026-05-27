using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Models: buoyancy_field

using ClimaSeaIce
using ClimaSeaIce: SeaIceModel, PhaseTransitions, ConductiveFlux
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium, PrescribedTemperature

using NumericalEarth.EarthSystemModels: ocean_surface_salinity
using NumericalEarth.Oceans
using NumericalEarth.SeaIces: sea_ice_dynamics


path_root="/data/hpcflash/users/josnez/Oceananigans/ICE-EDDY_wJ/V4/"


const Lx = 20kilometers
const Ly = 20kilometers
const Lz = 2kilometers

# Grid resolution
const Nx = 75
const Ny = 75
const Nz = 30

refinement = 10 
stretching = 10 

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

println("**************************")
println("Simulation set-up")

# =====================
# Grid (shared)
# =====================

using CUDA

Δt=5.0

grid = RectilinearGrid(GPU(),
                       size = (Nx, Ny, Nz),
                       x = (-Lx/2, Lx/2),
                       y = (-Ly/2, Ly/2),
                       halo = (7, 7, 7),
                       z = z_faces,
                       topology = (Bounded, Bounded, Bounded))


println("Ocean simulation")

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion=2.8e-4, haline_contraction=8e-4))

# Set up boundary conditions using Field
top_zonal_momentum_flux      = τˣ = Field{Face, Center, Nothing}(grid)
top_meridional_momentum_flux = τʸ = Field{Center, Face, Nothing}(grid)
top_ocean_heat_flux          = Jᵀ = Field{Center, Center, Nothing}(grid)
top_salt_flux                = Jˢ = Field{Center, Center, Nothing}(grid)

# Construct ocean boundary conditions including surface forcing and bottom drag
u_top_bc = FluxBoundaryCondition(τˣ)
v_top_bc = FluxBoundaryCondition(τʸ)
T_top_bc = FluxBoundaryCondition(Jᵀ)
S_top_bc = FluxBoundaryCondition(Jˢ)

boundary_conditions = (u = FieldBoundaryConditions(top=u_top_bc),
                       v = FieldBoundaryConditions(top=v_top_bc),
                       T = FieldBoundaryConditions(top=T_top_bc),
                       S = FieldBoundaryConditions(top=S_top_bc))


ocean_model = HydrostaticFreeSurfaceModel( grid;  buoyancy,
                            momentum_advection = WENO(),
                            tracer_advection = WENO(),
                            tracers = (:T, :S),
                            boundary_conditions=boundary_conditions)

#ocean = Simulation(ocean_model; Δt, verbose = false)

ocean = ocean_simulation(grid; Δt)

println("Ice simulation")

sea_ice = sea_ice_simulation(grid, ocean;
                            advection = WENO(order=7,
                            minimum_buffer_upwind_order=1),
                            )

#println("Atmosphere simulation")
#atmosphere = PrescribedAtmosphere(grid)


coupled_model = OceanSeaIceModel(ocean, sea_ice;
#                 ocean_reference_density=1026, 
#                 ocean_heat_capacity=3991.86795711963)
)

#coupled_model = EarthSystemModel(nothing, atmosphere, nothing, sea_ice, ocean;
#                                ocean_reference_density=1026,
#                                ocean_heat_capacity=3991.86795711963)

