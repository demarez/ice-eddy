using Oceananigans
using Oceananigans.Units
using CairoMakie
using Printf
using CUDA
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using NCDatasets
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
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
run_time = 181days
save_fields_interval = 24hour
path_root="/data/hpcflash/users/josnez/Oceananigans/ICE-EDDY_wJ/V0/"

###define paths
input_folder = path_root*"init_cond/"*experiment*"/" 
init_file = input_folder*"/init_julia.nc"
output_folder = path_root*"RUN/"*experiment*"/" 
mkpath(output_folder)

params = read_parameters_from_txt(input_folder*"parameters.txt")
println(params)


##### Dimensions of model: read in the parameters file from the init notebook

Lx = Int(params["Lx"]/1000)kilometers#200kilometers 
Ly = Int(params["Ly"]/1000)kilometers#200kilometers 
Lz = Int(params["Lz"]/1000)kilometers#2kilometers    

Nx = Int(params["ngrid_x"])#200
Ny = Int(params["ngrid_y"])#200
Nz = Int(params["ngrid_z"])#30

refinement = Int(params["refinement"])#4 
stretching = Int(params["stretching"])#40 

type_eddy=Int(params["type_eddy"])
cyclonic = type_eddy > 0

open("run.log", "w") do f
    # Write the parameters to the log file
    for (key, value) in params
        write(f, "$key $value\n")
    end
end

println("**************************")
println("GO")


rewrite = true#ARGS[2] #true #it has to be true the first time we run it !!!!
pickup = true#ARGS[3] #it has to be false if we want to re-write and re-start from beggining

##### ##### ##### ##### ##### ##### ##### 
##### Vertical spacing of Z ##### ##### 
##### ##### ##### ##### ##### ##### ##### 


# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)


##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Define grid, boundary and parameters ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### 

grid = RectilinearGrid(GPU(),
                       size = (Nx, Ny, Nz),
                       x = (-Lx/2, Lx/2),
                       y = (-Ly/2, Ly/2),
                       z = z_faces,
                       topology = (Bounded, Bounded, Bounded))

teos10 = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state=teos10)

#lineq_state = LinearEquationOfState(thermal_expansion=0, haline_contraction=8e-4)
#lineq_state = LinearEquationOfState(thermal_expansion=2.8e-4, haline_contraction=8.0e-4)
#buoyancy = SeawaterBuoyancy(equation_of_state=lineq_state)

no_slip_bc = ValueBoundaryCondition(0.0)
no_slip_field_bcs = FieldBoundaryConditions(no_slip_bc);

#-----#

cᴰ = 0 

if experiment == "idksomewithdrag"
    
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


# --- Edge Damping Setup ---
# A sponge layer is applied near the domain edges to damp velocities and prevent unphysical reflections.
# The damping is controlled by:
#   - EdgeMask{:xy}: Defines a smooth transition from full damping (A=1) at the edges to no damping (0) in the interior.
#     Parameters:
#       A=1: Amplitude of the mask.
#       f=0.255: Fraction of the domain width over which damping is applied.
#       delta=1/80: Sharpness of the transition between damped and undamped regions.
#       Lx, Ly: Domain extents in x and y directions.
#       threshold=0.001: Minimum mask value to avoid numerical issues.
#   - Relaxation: Applies a damping force proportional to the velocity, scaled by the mask.
#     rate=1/100: Relaxation timescale (100 time units).
# The damping is applied to u, v, and w velocities via the `forcing` keyword in the model.

A = 1
f = 0.255
delta = 1/80
threshold = 0.001
edge_mask = EdgeMask{:xy}(A=A, f=f, delta=delta, Lx=Lx, Ly=Ly, threshold=threshold )
damping = Relaxation(rate = 1/100, mask=edge_mask)
coriolis = FPlane(f=0.000143) #value at lat=80°N



if experiment == "idksomewithotherclosure"
    ν = 1e-4
    κ = 1e-4
    closure = ScalarBiharmonicDiffusivity(; ν, κ)
    #closure = ScalarDiffusivity(; ν, κ)
    println("ScalarBiharmonicDiffusivity")
#elseif experiment == "test1_window"
    # Do something else
#    println("Running test1_window")
else
    closure = CATKEVerticalDiffusivity()
    println("CATKEVerticalDiffusivity")
end




model = HydrostaticFreeSurfaceModel(; grid, buoyancy,
                            momentum_advection = WENO(),
                            tracer_advection = WENO(),
                            tracers = (:T, :S, :e),
                            coriolis = coriolis,
                            closure = closure,
                            boundary_conditions=(u=u_bcs, v=v_bcs),
                            forcing=(u=damping, v=damping, w=damping)
)

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

# `set!` the `model` fields using functions or constants:
set!(model, u=uᵢ, v=vᵢ, T=Tᵢ, S=Sᵢ, η=ηᵢ)


max_dt = 30second # using automatic submit, this should be in seconds

simulation = Simulation(model, Δt=1.0, stop_time=run_time)

wizard = TimeStepWizard(cfl=1.0, max_change=1.1, max_Δt=max_dt)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

function progress(simulation)
    u, v, w = simulation.model.velocities

    # Print a progress message
    msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, wall time: %s\n",
                   iteration(simulation),
                   prettytime(time(simulation)),
                   prettytime(simulation.Δt),
                   maximum(u), maximum(v), maximum(w),
                   prettytime(simulation.run_wall_time))

    @info msg
    
    # Write to run.log file (append mode)
    open("run.log", "a") do f
        write(f, msg)
    end

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))

u, v, w = model.velocities

T,S = model.tracers
b = BuoyancyField(model)

ζ = ∂x(v) - ∂y(u)



tracer_fields = Dict("T" => T, "S" => S, "B" => b);

vel_fields = Dict("U" => u, "V" => v, "vort" => ζ);

println("create output fields")

simulation.output_writers[:tracer_field_writer] =
    NetCDFOutputWriter(model, tracer_fields; 
                       filename=output_folder*"tracer_fields.nc", 
                       schedule = TimeInterval(save_fields_interval),
                       overwrite_existing = rewrite,
#                       file_splitting = TimeInterval(30days),
                       )

simulation.output_writers[:vel_field_writer] =
    NetCDFOutputWriter(model, vel_fields; 
                       filename=output_folder*"vel_fields.nc", 
                       schedule = TimeInterval(save_fields_interval),
                       overwrite_existing = rewrite,
#                       file_splitting = TimeInterval(30days),
                       )

simulation.output_writers[:checkpointer] = 
    Checkpointer(model;
		 dir = output_folder,
                 schedule=TimeInterval(30days), 
                 prefix="model_checkpoint",
                 overwrite_existing = rewrite,
                 verbose = true
                 )
println("RUN!")
run!(simulation, pickup=pickup)
