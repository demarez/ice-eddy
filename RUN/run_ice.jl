# Idealized Ocean-Ice coupled simulation.
# This script buids from the PR #54 of Simone-Silvestri: Wind-Driven Ice-Ocean Channel: AB2+FE vs RK3+RK3 Coupling Comparison
# ===========================================================================

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Models: buoyancy_field

using ClimaSeaIce
using ClimaSeaIce: SeaIceModel, SlabSeaIceThermodynamics, PhaseTransitions, ConductiveFlux
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium, PrescribedTemperature

using NumericalEarth.EarthSystemModels: ocean_surface_salinity
using NumericalEarth.Oceans

using NCDatasets
using Printf
using Base.Filesystem

include("./utils.jl")


#a function to transform string to bool
parsebool(s::String) = lowercase(s) == "true" ? true : false


function read_parameters_from_txt(path_read::String)
    params = Dict{String, Float64}()

    open(path_read, "r") do file
        for line in eachline(file)
            parts = split(line)
            if length(parts) >= 2
                param_name = parts[1]
                param_value = tryparse(Float64, parts[2])
                if param_value !== nothing
                    params[param_name] = param_value
                else
                    # Handle non-numeric parameters (e.g., type_eddy)
                    params[param_name] = join(parts[2:end], " ")
                end
            end
        end
    end

    return params
end

experiment = ARGS[1]
run_time = 45days#91days#181days
save_fields_interval = 24hour
path_root="/data/hpcflash/users/josnez/Oceananigans/ICE-EDDY_wJ/V4/"

###define paths
input_folder = path_root*"init_cond/"*experiment*"/" 
init_file = input_folder*"/init_julia.nc"
output_folder = path_root*"RUN/"*experiment*"/" 
mkpath(output_folder)

params = read_parameters_from_txt(input_folder*"parameters.txt")
println(params)


# =====================
# Physical parameters
# =====================

# Domain

const Lx = Float32(params["Lx"]/1000)kilometers#200kilometers 
const Ly = Float32(params["Ly"]/1000)kilometers#200kilometers 
const Lz = Float32(params["Lz"]/1000)kilometers#2kilometers  

# Grid resolution
const Nx = Int(params["ngrid_x"])#200
const Ny = Int(params["ngrid_y"])#200
const Nz = Int(params["ngrid_z"])#30

refinement = Int(params["refinement"])#4 
stretching = Int(params["stretching"])#40 

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

type_eddy=Int(params["type_eddy"])
cyclonic = type_eddy > 0

open("run_"*experiment*".log", "w") do f
    # Write the parameters to the log file
    for (key, value) in params
        write(f, "$key $value\n")
    end
end

println("**************************")
println("GO")

rewrite = true#ARGS[2] #true #it has to be true the first time we run it !!!!
pickup = true#ARGS[3] #it has to be false if we want to re-write and re-start from beggining

# =====================
# Grid (shared)
# =====================

using CUDA

grid = RectilinearGrid(GPU(),
                       size = (Nx, Ny, Nz),
                       x = (-Lx/2, Lx/2),
                       y = (-Ly/2, Ly/2),
                       halo = (7, 7, 7),
                       z = z_faces,
                       topology = (Bounded, Bounded, Bounded))

# =====================
# Damping region
# =====================

A = 1
f = 0.255
delta = 1/80
threshold = 0.001
edge_mask = EdgeMask{:xy}(A=A, f=f, delta=delta, Lx=Lx, Ly=Ly, threshold=threshold )
damping = Relaxation(rate = 1/1000, mask=edge_mask)

# =====================
# Initial conditions
# =====================

# Sea ice initial conditions
const h₀ = 0.5   # initial ice thickness (m)
const ℵ₀ = 0.9   # initial ice concentration

# Ocean initial conditions

ds = Dataset(init_file);

get_data(var) = Float64.(ds[var][:,:,end:-1:1])
get_ssh(var) = Float64.(ds[var][:,:])

# Temperature and salinity initial condition: a stable density gradient with random noise superposed.
Tᵢ = get_data("votemper")
Sᵢ = get_data("vosaline")

ηᵢ = get_ssh("ssh")

# Velocity initial condition.
uᵢ = get_data("uoce")
vᵢ = get_data("voce")

# =====================
# Model parameters
# =====================

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion=2.8e-4, haline_contraction=8e-4))

# Coriolis (This is northern hemisphere)
const f₀ = -0.000143  # s⁻¹
coriolis = FPlane(f=f₀)

ν = 1e-4
κ = 1e-4
closure = ScalarBiharmonicDiffusivity(; ν, κ)
println("ScalarBiharmonicDiffusivity")

# =====================
# Model simulation
# =====================

output_prefix = "OceanIce_periodic"

Δt=5.0
max_dt = 30second # using automatic submit, this should be in seconds
run_time = 91days

# --- Atmosphere (fresh instance) ---
# Atmospheric state (polar winter)
const Tₐ  = 273.15 - 20  # -20°C in Kelvin
const u₁₀ = 0            # No winds
const qₐ  = 0.001        # specific humidity (dry cold air)
const Qsw = 0.0          # no shortwave (polar winter)
const Qlw = 180.0        # downwelling longwave (W/m²)

function build_atmosphere(arch)
    atmosphere_grid  = RectilinearGrid(arch; size = (1, 1),
                                       x = (0, Lx),
                                       y = (0, Ly),
                                       topology = (Periodic, Bounded, Flat))
    atmosphere_times = range(0, 365days, length=3)
    atmosphere = PrescribedAtmosphere(atmosphere_grid, atmosphere_times)

    parent(atmosphere.tracers.T)  .= Tₐ
    parent(atmosphere.velocities.u) .= u₁₀
    parent(atmosphere.tracers.q)  .= qₐ
    # parent(atmosphere.downwelling_radiation.shortwave) .= Qsw
    # parent(atmosphere.downwelling_radiation.longwave)  .= Qlw

    return atmosphere
end

atmosphere = build_atmosphere(Oceananigans.Architectures.architecture(grid))
radiation  = Radiation(ocean_albedo=0.06, sea_ice_albedo=0.7)

cᴰ = 0

if occursin("drag",experiment)
    z₀ = 0.01 # m (roughness length) ###the one we vary
    κ = 0.4  # von Karman constant
    z₁ = CUDA.@allowscalar -last(znodes(grid, Center())) # Closest grid center to the bottom
    # Drag coefficient
    cᴰ = (κ / log(z₁ / z₀))^2
    println("some drag")
#elseif experiment == "test1_window"
    # Do something else
#    println("Running test1_window")
else
    cᴰ = 0
    println("No drag")
end

Uᵢ = 0 # m s⁻¹
Vᵢ = 0 # m s⁻¹
rho = 1026.0

@inline drag_u(x, y, t, u, v, p) =  p.rho * p.cᴰ * √((u - p.Uᵢ)^2 + (v - p.Vᵢ)^2) * ( u - p.Uᵢ)
@inline drag_v(x, y, t, u, v, p) =  p.rho * p.cᴰ * √((u - p.Uᵢ)^2 + (v - p.Vᵢ)^2) * ( v - p.Vᵢ)

drag_bc_u = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v), parameters=(; cᴰ, Uᵢ, Vᵢ, rho))
drag_bc_v = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v), parameters=(; cᴰ, Uᵢ, Vᵢ, rho))

u_bcs = FieldBoundaryConditions( top = drag_bc_u )
v_bcs = FieldBoundaryConditions( top = drag_bc_v )


# --- Ocean ---
#ocean = ocean_simulation(grid;
#                            Δt,
#                            coriolis = coriolis,
#                            closure = closure,
#                            momentum_advection = WENOVectorInvariant(),
#                            tracer_advection = WENO(order=7),
#                            forcing=(u=damping, v=damping))

bottom_drag_coefficient = 0.003

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
                            coriolis = coriolis,
                            closure = closure,
                            boundary_conditions=boundary_conditions,
                            forcing=(u=damping, v=damping))

ocean = Simulation(ocean_model; Δt, verbose = false)

#wizard = TimeStepWizard(cfl=0.5, max_change=1.1, max_Δt=max_dt)
#ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

set!(ocean.model, u=uᵢ, v=vᵢ, T=Tᵢ, S=Sᵢ, η=ηᵢ)

# sea_ice_dyn = sea_ice_dynamics(grid, ocean; sea_ice_ocean_drag_coefficient = 3.24e-3)

sea_ice = sea_ice_simulation(grid, ocean; 
                            advection = WENO(order=7, 
                            minimum_buffer_upwind_order=1))#,
                            # dynamics = sea_ice_dyn)
                            
set!(sea_ice.model, h=h₀, ℵ=ℵ₀)

println("Set up coupled model")

const ocean_rho_ref = 1030
const ocean_heat_cap = 3991.86795711963

# --- Coupled Model ---
coupled_model = OceanSeaIceModel(ocean, sea_ice; 
                                 atmosphere, radiation, 
                                 ocean_reference_density=ocean_rho_ref, 
                                 ocean_heat_capacity=ocean_heat_cap)

#coupled_model = OceanSeaIceModel(ocean, sea_ice)

simulation    = Simulation(coupled_model; Δt, stop_time=run_time)

function progress(simulation)
    u, v, w = ocean.model.velocities

    # Print a progress message
    msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, wall time: %s\n",
                   iteration(simulation),
                   prettytime(time(simulation)),
                   prettytime(simulation.Δt),
                   maximum(u), maximum(v), maximum(w),
                   prettytime(simulation.run_wall_time))

    @info msg
    
    # Write to run.log file (append mode)
    open("run_"*experiment*".log", "a") do f
        write(f, msg)
    end

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))
u, v, w = ocean.model.velocities

T,S = ocean.model.tracers
b = buoyancy_field(ocean.model)

ζ = ∂x(v) - ∂y(u)

tracer_fields = Dict("T" => T, "S" => S, "B" => b);

vel_fields = Dict("U" => u, "V" => v, "vort" => ζ);

println("create output fields")

simulation.output_writers[:tracer_field_writer] =
             NetCDFWriter(ocean.model, tracer_fields; 
                       filename=output_folder*"tracer_fields.nc", 
                       schedule = TimeInterval(save_fields_interval),
                       overwrite_existing = rewrite,
#                       file_splitting = TimeInterval(30days),
                       )

simulation.output_writers[:vel_field_writer] =
             NetCDFWriter(ocean.model, vel_fields; 
                       filename=output_folder*"vel_fields.nc", 
                       schedule = TimeInterval(save_fields_interval),
                       overwrite_existing = rewrite,
#                       file_splitting = TimeInterval(30days),
                       )

sea_ice_outputs = (; h = sea_ice.model.ice_thickness,
                     ℵ = sea_ice.model.ice_concentration)

simulation.output_writers[:ice_fields_writer] = 
             NetCDFWriter(sea_ice.model, sea_ice_outputs;
                       filename=output_folder*"vel_fields.nc",
                       schedule = TimeInterval(save_fields_interval),
                       overwrite_existing = rewrite,
                       )

# Ocean net surface fluxes (momentum + tracers)
ocean_flux_outputs = (; τx = coupled_model.interfaces.net_fluxes.ocean.u,
                        τy = coupled_model.interfaces.net_fluxes.ocean.v,
                        JT = coupled_model.interfaces.net_fluxes.ocean.T,
                        JS = coupled_model.interfaces.net_fluxes.ocean.S)

simulation.output_writers[:fluxes] = 
             NetCDFWriter(ocean.model, ocean_flux_outputs;
                       filename = output_folder*"ocean_fluxes.nc",
                       schedule = TimeInterval(save_fields_interval),
                       overwrite_existing = rewrite,
                       )

simulation.output_writers[:checkpointer] = 
    Checkpointer(ocean.model;
		 dir = output_folder,
                 schedule=TimeInterval(30days), 
                 prefix="model_checkpoint",
                 overwrite_existing = rewrite,
                 verbose = true
                 )

println("RUN!")
run!(simulation)#, pickup=pickup)
