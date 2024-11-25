# fluid-sim-1
this represents my attempt at creating a fluid mechanics discrete solver
it has the following models:
	incompressible momentum solving,
	heat diffusion,
	heat convection,
	buoyancy modeled using the Boussinesq approximation
each model is enabled/disabled from the init file (an example is included)
solver parameters (time step size, length scale, etc.), as well as fluid and solid properties are also specified in the init file

the syntax of each init file item are explained below

SOLVER PARAMETERS
	WIDTH [int] - the width, in cells, of the domain
	HEIGHT [int] - the height, in cells, of the domain
	L_SCALE [float] - the width, in arbitrary length units, of each cell
	T_STEP [float or fraction] - the length, in seconds, of each time step (a fraction is formatted as [x/y])
	ANIM_FREQ [int] - the number of time steps completed before rendering a frame
	SAVE_FREQ [int] - the number of frames rendered before saving data to file
	ERR_THRESH [float] - the total amount of momentum error allowed before progressing to the next time step
	MAX_ITERS [int] - the maximum number of iterations performed for each time step before progressing (if residuals are still above the allowed maximum)
	G [float] - the gravitational constant (used only if buoyancy is enabled)

MODEL ENABLE (all have [1] or [0] as an argument for enabled or disabled, respectively. additional arguments are list below)
	THERM_DIFFUSE [int] - the number of time subdivisions to perform heat equations on each time step
	MOMENTUM RELAX [float] - if enabled, a relax constant (rec between 1.0 and 1.9) is also specified
	CONVECTION
	BUOYANCY
	VISCOUS

PROPERTIES
	RHO is the density, T0 the starting temperature, P0 the starting pressure,
	KAPPA the conductivity, ALPHA the coefficient of thermal expansion,
	C the heat capactiy, and MU the viscocity
	fluid and solid phases have their accepted arguments listed below

	FLUID RHO [float] T0 [float] KAPPA [float] ALPHA [float] C [float] MU [float]
	SOLID RHO [float] KAPPA [float] C [float]

BOUNDARY CONDITIONS (currently, only wall boundary conditions are implemented)
	[TOP/BOTTOM/RIGHT/LEFT] [WALL] [AD/TEMP/FLUX] [VALUE]
	walls can be AD (adiabatic), TEMP (a fixed temperature), or FLUX (a heat flux)
	if not adiabatic, the value of the temperature or flux must be specified

ADDITIONAL FLUXES
	FLUX [TOP/BOTTOM/LEFT/RIGHT] [int:starting square] [int:ending square] [value]

INITIAL CONDITONS
	COND [TEMP] [float:value] [RECT/CIRCLE] [list:dimensions]
	for a rectangle, the dimensions are [int:corner x] [int:corner y] [int:width] [int:height]
	for a circle, the dimensions are [int:center x] [int:center y] [int: radius]
	all dimensions are in cells, not in absolute units
