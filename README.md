# 2D fluids solver/sandbox using pillow library for vizualization

import sys, math, os
from PIL import Image, ImageDraw, ImageFont

# a few basic color values
AQUA = ( 0, 255, 255)
BLACK = ( 0, 0, 0)
BLUE = ( 0, 0, 255)
FUCHSIA = (255, 0, 255)
GREY = (128, 128, 128)
GREEN = ( 0, 128, 0)
LIME = ( 0, 255, 0)
MAROON = (128, 0, 0)
NAVY = ( 0, 0, 128)
OLIVE = (128, 128, 0)
PURPLE = (128, 0, 128)
RED = (255, 0, 0)
SILVER = (210, 210, 210)
TEAL = ( 0, 128, 128)
WHITE = (255, 255, 255)
YELLOW = (255, 255, 0)

MAX_FRAMES = 10000

MAX_ITERS = 2

def main():
    global PARAMS
    global MODELS
    global EDGES
    global PROPERTIES
    global FLUXES
    global PORTS
    global CONDS
    PARAMS, MODELS, EDGES, PROPERTIES, FLUXES, PORTS, CONDS = init_params()

    saveFile = input("Enter save file name to load data, or press enter to re-initialize")

    global solids

    if saveFile == "":
        temps = []
        pressures = []
        solids = []

        # IMPORTANT: the j coordinate (ie, the column) will be used for the points that can interact with each other
        # this means that the horizontal edge points will be column, row (ie, horts[col][row])
        verts = []
        horts = []

        for i in range(PARAMS["HEIGHT"]):
            verts.append([])
            for j in range(PARAMS["WIDTH"] + 1):
                verts[i].append(0)

        for i in range(PARAMS["WIDTH"]):
            horts.append([])
            for j in range(PARAMS["HEIGHT"] + 1):
                horts[i].append(0)

        # generate grid squares
        for i in range(PARAMS["HEIGHT"]):
            temps.append([])
            pressures.append([])
            solids.append([])
            for j in range(PARAMS["WIDTH"]):
                temps[i].append(PROPERTIES["FLUID"]["T0"])
                pressures[i].append(PROPERTIES["FLUID"]["P0"])
                solids[i].append(1)

        for cond in CONDS:
            changeSquares = []
            if cond["SHAPE"] == "CIRCLE":
                centerX = int(cond["MEASURE"][0])
                centerY = int(cond["MEASURE"][1])
                radius = int(cond["MEASURE"][2])

                for i in range(PARAMS["HEIGHT"]):
                    for j in range(PARAMS["WIDTH"]):
                        if math.sqrt((i - centerY)**2 + (j - centerX)**2) < radius:
                            changeSquares.append([i,j])

            elif cond["SHAPE"] == "RECT":
                cornerX = int(cond["MEASURE"][0])
                cornerY = int(cond["MEASURE"][1])
                width = int(cond["MEASURE"][2])
                height = int(cond["MEASURE"][3])

                for i in range(PARAMS["HEIGHT"]):
                    for j in range(PARAMS["WIDTH"]):
                        if i > cornerY and i < cornerY + height:
                            if j > cornerX and j < cornerX + width:
                                changeSquares.append([i,j])
                
            if cond["TYPE"] == "TEMP":
                for sq in changeSquares:
                    temps[sq[0]][sq[1]] = cond["QUANT"]

        t_elapsed = 0

        frames = 0

        resid = 0

        minT, maxT = getBounds(temps)
        minP, maxP = getBounds(pressures)

        generateImage(temps, minT, maxT, pressures, minP, maxP, resid, t_elapsed, frames)

    else:
        newData = read_data(saveFile)

        for label in newData:
            print(label)

        temps = newData["temps"]
        pressures = newData["pressures"]

        
        solids = []

        for i in range(PARAMS["HEIGHT"]):
            solids.append([])
            for j in range(PARAMS["WIDTH"]):
                solids[i].append(True)
                
        verts = newData["verts"]
        horts = newData["horts"]

        t_elapsed = newData["t_elapsed"][0]
        frames = int(newData["frames"][0])

        print("Resuming from frame", frames, "t =", t_elapsed, "s")

    keep_going = True

    sinceLastSave = 1

    #fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", 40)


    # this array will be used for saving processor resources by not updating every square at each time step
    # any square that gets updated will put all its neighbors as needing to be updated on the next iteration
    # also any square with buoyancy effects applied will also be updated in momentum
    dont_calc = []

    for i in range(PARAMS["HEIGHT"]):
            dont_calc.append([])
            for j in range(PARAMS["WIDTH"] + 1):
                dont_calc[i].append(False)
    
    # <--------------------------------- main game loop -------------------->
    while keep_going:

        for frame in range(PARAMS["ANIM_FREQ"]):
            if MODELS["THERM_DIFFUSE"]:
                temps = heat_diffuse(temps)
            if MODELS["VISCOUS"]:
                horts, verts = apply_viscocity(horts, verts)
            if MODELS["BUOYANCY"]:
                horts, dont_calc = apply_buoyancy(horts, temps, dont_calc)
            if MODELS["MOMENTUM"]:
                iters = 1
                horts, verts, pressures, resid, dont_calc = solve_continuity(horts, verts, pressures, dont_calc)

                while resid > PARAMS["ERR_THRESH"] and iters < PARAMS["MAX_ITERS"]:
                    horts, verts, pressures, resid, dont_calc = solve_continuity(horts, verts, pressures, dont_calc)
                    iters += 1

                print(resid, iters, frame)
                
            if MODELS["CONVECTION"]:
                temps = solve_convection(horts, verts, temps)

            t_elapsed += PARAMS["T_STEP"]

        frames += 1

        print("frame",frames,"t elapsed:",t_elapsed,"s")

        minT, maxT = getBounds(temps)
        minP, maxP = getBounds(pressures)

        generateImage(temps, minT, maxT, pressures, minP, maxP, resid, t_elapsed, frames)

        if sinceLastSave < PARAMS["SAVE_FREQ"]:
            sinceLastSave += 1

        else:
            save_data("save"+str(frames)+".txt", ["temps", "pressures", "verts", "horts", "t_elapsed", "frames"],
                      [2, 2, 2, 2, 1, 1], [temps, pressures, verts, horts, [t_elapsed,], [frames,] ] )

            print("frame", frames, "saved, t=", t_elapsed, "s")
            sinceLastSave = 1         

        if frames >= MAX_FRAMES:
            keep_going = False

def pres_to_color(p, minP, maxP):
    if minP != maxP:
        relP = (p - minP) / (maxP - minP)
    else:
        relP = 0.5
    if p >= maxP:
        return (255,255,255)
    
    a = int(-255 * math.sqrt(1 - relP) + 255)

    if a < 0:
        return (0,0,0)

    if a > 255:
        return (255,255,255)

    return (a,a,a)

def temp_to_color(t, minT, maxT):
    a = -1/255
    if maxT != minT:
        tRel = ( t - minT) / (maxT - minT)
    else:
        tRel = 0.5
    r = int(-255 * math.sqrt(1 - tRel) + 255)
    try:
        g = int(-4 * 255 * tRel * (tRel - 1))
    except:
        g = 255
    try:
        b = int(255 * (1 - math.sqrt(tRel)))
    except:
        b = 255

    if r < 0:
        r = 0
    elif r > 255:
        r = 255

    if g < 0:
        g = 0
    elif g > 255:
        g = 255

    if b < 0:
        b = 0
    elif b > 255:
        b = 255

    return (r,g,b)

def heat_diffuse(temps):
    newTemps = temps.copy()
    oldTemps = newTemps.copy()

    for a in range(PARAMS["TEMP_ITERS"]):
        const = PARAMS["T_STEP"] / (PROPERTIES["FLUID"]["C"] * PROPERTIES["FLUID"]["RHO"] * PARAMS["L_SCALE"] * PARAMS["TEMP_ITERS"])

        # handle other fluxes
        for newFlux in FLUXES:
            if newFlux["EDGE"] == "TOP":
                for j in range(newFlux["START"], newFlux["END"] + 1):
                    newTemps[0][j] += newFlux["QUANT"] * const

            elif newFlux["EDGE"] == "BOTTOM":
                for j in range(newFlux["START"], newFlux["END"] + 1):
                    newTemps[-1][j] += newFlux["QUANT"] * const

            elif newFlux["EDGE"] == "LEFT":
                for i in range(newFlux["START"], newFlux["END"] + 1):
                    newTemps[i][0] += newFlux["QUANT"] * const

            elif newFlux["EDGE"] == "RIGHT":
                for i in range(newFlux["START"], newFlux["END"] + 1):
                    newTemps[i][-1] += newFlux["QUANT"] * const
        
        for i in range(PARAMS["HEIGHT"]):
            for j in range(PARAMS["WIDTH"]):
                if solids[i][j] == 0:
                    k = PROPERTIES["SOLID"]["KAPPA"]
                    hc = PROPERTIES["SOLID"]["C"]
                    density = PROPERTIES["SOLID"]["RHO"]

                else:
                    k = PROPERTIES["FLUID"]["KAPPA"]
                    hc = PROPERTIES["FLUID"]["C"]
                    density = PROPERTIES["FLUID"]["RHO"]

                # law of heat diffusion is q = -k dT / dx
                # the resulting change in temperature is the heat flux times square width times time step size divided by the specific heat capacity times the area times the density
                # length scale should cancel out because the heat flux per unit length perpendicular to the plane
                # is the heat flux per unit area multiplied by the length scale
                # this cancels out the length scale in the denominator
                denom = hc * PARAMS["L_SCALE"] ** 2 * density * PARAMS["TEMP_ITERS"]

                const = k * PARAMS["T_STEP"] / denom

                # top edge
                if i == 0:
                    if EDGES["TOP"]["TYPE"] == "WALL":
                        if EDGES["TOP"]["AD"]:
                            newTemps[i][j] += const * (oldTemps[i+1][j] - oldTemps[i][j]) # adiabatic wall, only the non-edge-square conducts heat

                        elif EDGES["TOP"]["TEMP"]:
                            newTemps[i][j] += const * (EDGES["TOP"]["QUANT"] + oldTemps[i+1][j] - 2 * oldTemps[i][j]) # temperature boundary, sub in the boundary temp for a square
                            
                        elif EDGES["TOP"]["FLUX"]:
                            newTemps[i][j] += const * (oldTemps[i+1][j] - oldTemps[i][j] + EDGES["TOP"]["QUANT"] * PARAMS["L_SCALE"]) # heat flux boundary, conducts with non-edge, add the heat flux (divided by k because there's a k in the numerator)
                        
                # right edge
                elif i == PARAMS["HEIGHT"] - 1:
                    if EDGES["BOTTOM"]["TYPE"] == "WALL":
                        if EDGES["BOTTOM"]["AD"]:
                            newTemps[i][j] += const * (oldTemps[i-1][j] - oldTemps[i][j])

                        elif EDGES["BOTTOM"]["TEMP"]:
                            newTemps[i][j] += const * (EDGES["BOTTOM"]["QUANT"] + oldTemps[i-1][j] - 2 * oldTemps[i][j])

                        elif EDGES["BOTTOM"]["FLUX"]:
                            newTemps[i][j] += const * (oldTemps[i-1][j] - oldTemps[i][j] + EDGES["BOTTOM"]["QUANT"] * PARAMS["L_SCALE"])

                else:
                    newTemps[i][j] += const * (oldTemps[i-1][j] + oldTemps[i+1][j] - 2 * oldTemps[i][j])

                # left edge
                if j == 0:
                    if EDGES["LEFT"]["TYPE"] == "WALL":
                        if EDGES["LEFT"]["AD"]:
                            newTemps[i][j] += const * (oldTemps [i][j+1] - oldTemps[i][j])

                        elif EDGES["LEFT"]["TEMP"]:
                            newTemps[i][j] += const * (EDGES["LEFT"]["QUANT"] + oldTemps [i][j+1] - 2 * oldTemps[i][j])

                        elif EDGES["LEFT"]["FLUX"]:
                            newTemps[i][j] += const * (oldTemps [i][j+1] - oldTemps[i][j] + EDGES["LEFT"]["QUANT"] * PARAMS["L_SCALE"])
                        

                # right edge
                elif j == PARAMS["WIDTH"] - 1:
                    if EDGES["RIGHT"]["TYPE"] == "WALL":
                        if EDGES["RIGHT"]["AD"]:
                            newTemps[i][j] += const * (oldTemps[i][j-1] - oldTemps[i][j])

                    elif EDGES["RIGHT"]["TEMP"]:
                        newTemps[i][j] += const * (EDGES["RIGHT"]["QUANT"] + oldTemps[i][j-1] - 2 * oldTemps[i][j])
                        
                    elif EDGES["RIGHT"]["FLUX"]:
                        newTemps[i][j] += const * (oldTemps[i][j-1] - oldTemps[i][j] + EDGES["RIGHT"]["QUANT"] * PARAMS["L_SCALE"])

                else:
                    newTemps[i][j] += const * (oldTemps[i][j-1] + oldTemps[i][j+1] - 2 * oldTemps[i][j])

        oldTemps = newTemps.copy()

    return newTemps

# function to solve continuity equations
# the sum of mass flows into a square should be zero (divergence free flow field)
def solve_continuity(horts, verts, pressures, dont_calc):
    newHorts = horts.copy()
    newVerts = verts.copy()
    newPressures = pressures.copy()

    new_dont_calc = []

    for i in range(PARAMS["HEIGHT"]):
            new_dont_calc.append([])
            for j in range(PARAMS["WIDTH"] + 1):
                new_dont_calc[i].append(True)

    #print("Going around again")
    total_err = 0
    for i in range(PARAMS["HEIGHT"]):
        for j in range(PARAMS["WIDTH"]):
            # dont calculate for this square if deemed unnecessary
            if dont_calc[i][j]:
                continue
            
            if solids[i][j] == 0:
                continue
            
            momentum_err = 0
            s = 0 # number of non-controlled edges of the square

            # top edge
            if i == 0:
                if EDGES["TOP"]["TYPE"] == "WALL":
                    s += 0
                # top BC velocity minus velocity at bottom of square
                #if TOP_V != "copy":
                #    momentum_err += TOP_V - horts[j][i+1]
                #    s += 2

            # not top edge
            else:
                s += solids[i-1][j]
                momentum_err += horts[j][i]

            # bottom edge
            if i == PARAMS["HEIGHT"] - 1:
                if EDGES["BOTTOM"]["TYPE"] == "WALL":
                    s += 0
                # velocity at top of square minus bottom BC velocity
                #if BOT_V != "copy":
                #    momentum_err += horts[j][i] - BOT_V
                #    s += 2

            # not bottom edge
            else:
                s += solids[i+1][j]
                momentum_err += -horts[j][i+1]

            # left edge
            if j == 0:
                if EDGES["LEFT"]["TYPE"] == "WALL":
                    s += 0
                # left BC velocity minus velocity at right of square
                #if LEF_V != "copy":
                #    momentum_err += LEF_V - verts[i][j+1]
                #    s += 2

            # not on left edge
            else:
                s += solids[i][j-1]
                momentum_err += verts[i][j]

            # right edge
            if j == PARAMS["WIDTH"] - 1:
                if EDGES["RIGHT"]["TYPE"] == "WALL":
                    s += 0
                # velocity at left of square minus right BC velocity
                #if RIG_V != "copy":
                #    momentum_err += verts[i][j] - RIG_V
                #    s += 2
                    
            # not on right edge
            else:
                s += solids[i][j+1]
                momentum_err += -verts[i][j+1]
                        
            momentum_adj = momentum_err * PARAMS["RELAX"] / s

            # if the momentum within square i,j had to be adjusted, make sure that it and all its neighbors are adjusted on the next iteration
            if abs(momentum_adj) > 0:
                new_dont_calc[i][j] = False

                if i > 0:
                    new_dont_calc[i-1][j] = False
                if i < PARAMS["HEIGHT"] - 1:
                    new_dont_calc[i+1][j] = False
                if j > 0:
                    new_dont_calc[i][j-1] = False
                if j < PARAMS["WIDTH"] - 1:
                    new_dont_calc[i][j+1] = False
            
            # adjust each edge of the square
            # but not domain edges
            if i > 0:
                newHorts[j][i] += -momentum_adj * solids[i-1][j] # top
            if i < PARAMS["HEIGHT"]-1:
                newHorts[j][i+1] += momentum_adj * solids[i+1][j] # bottom
            if j > 0:
                newVerts[i][j] += -momentum_adj * solids[i][j-1] # left
            if j < PARAMS["WIDTH"]-1:
                newVerts[i][j+1] += momentum_adj * solids[i][j+1] # right

            newPressures[i][j] += momentum_adj * PROPERTIES["FLUID"]["RHO"] * PARAMS["L_SCALE"] / (s * PARAMS["T_STEP"])
            
            total_err += abs(momentum_err * PROPERTIES["FLUID"]["RHO"] * PARAMS["L_SCALE"])

    # then, re-set boundary conditions
    for i in range(len(verts)):
        if EDGES["LEFT"]["TYPE"] == "WALL":
            newVerts[i][0] = 0

        if EDGES["RIGHT"]["TYPE"] == "WALL":
            newVerts[i][-1] = 0

    for i in range(len(horts)):
        if EDGES["TOP"]["TYPE"] == "WALL":
            newHorts[i][0] = 0

        if EDGES["BOTTOM"]["TYPE"] == "WALL":
            newHorts[i][-1] = 0

    dom_mass = PARAMS["HEIGHT"] * PARAMS["WIDTH"] * PARAMS["L_SCALE"] ** 2 * PROPERTIES["FLUID"]["RHO"]

    resid = total_err / dom_mass

    return newHorts, newVerts, newPressures, resid, new_dont_calc

def apply_buoyancy(horts, temps, dont_calc):
    newHorts = horts.copy()
    new_dont_calc = dont_calc.copy()

    const = 0.5 * PROPERTIES["FLUID"]["ALPHA"] * PARAMS["G"] * PARAMS["T_STEP"] / PROPERTIES["FLUID"]["RHO"]

    # apply effects of buoyancy
    for j in range(PARAMS["WIDTH"]):
        if temps[0][j] != PROPERTIES["FLUID"]["T0"]:
            new_dont_calc[0][j] = False
            newHorts[j][1] += const * (PROPERTIES["FLUID"]["T0"] - temps[0][j]) # apply effects of top cell of each column to the value second down
        if temps[-1][j] != PROPERTIES["FLUID"]["T0"]:
            new_dont_calc[-1][j] = False
            newHorts[j][-1] += const * (PROPERTIES["FLUID"]["T0"] - temps[-1][j]) # apply effects of the bottom cell to the value second up

        # go thru the rest of the column
        # because each velocity value will be added to twice, multiply the effect of each cell by half
        for i in range(1,PARAMS["HEIGHT"]-1):
            dvx = const * (PROPERTIES["FLUID"]["T0"] - temps[i][j])
            if dvx != 0:
                new_dont_calc[i][j] = False
                newHorts[j][i] += dvx
                newHorts[j][i+1] += dvx

    return newHorts, new_dont_calc

def apply_viscocity(horts, verts):
    newHorts = horts.copy()
    newVerts = verts.copy()

    denom = PARAMS["L_SCALE"] ** 3 * PROPERTIES["FLUID"]["RHO"]

    const = PROPERTIES["FLUID"]["MU"] * PARAMS["T_STEP"] / denom

    for i in range(PARAMS["HEIGHT"]):
        for j in range(PARAMS["WIDTH"]):
            if solids[i][j] == 0:
                continue
            
            u_net = verts[i][j+1] - verts[i][j] # net x-velocity thru square
            v_net = horts[j][i-1] - horts[j][i] # net y-velocity thru square
            
            if i == 0:
                if EDGES["TOP"]["TYPE"] == "WALL":
                    du_up = 0
                else:
                    du_up = u_net

            else:
                du_up = verts[i-1][j+1] - verts[i-1][j] # net x-vel thru square above

            if i == PARAMS["HEIGHT"] - 1:
                if EDGES["BOTTOM"]["TYPE"] == "WALL":
                    du_down = 0
                else:
                    du_down = u_net

            else:
                du_down = verts[i+1][j+1] - verts[i+1][j] # net thru square below

            newVerts[i][j] += const * (du_down + du_up - 2 * u_net)
            newVerts[i][j+1] += const * (du_down + du_up * 2 * u_net)

            if j == 0:
                if EDGES["LEFT"]["TYPE"] == "WALL":
                    dv_left = 0
                else:
                    dv_left = v_net

            else:
                dv_left = horts[j][i+1] - horts[j][i+1]

            if j == PARAMS["WIDTH"] - 1:
                if EDGES["RIGHT"]["TYPE"] == "WALL":
                    dv_right = 0
                else:
                    dv_right = v_net

            else:
                dv_right = horts[j+1][i+1] - horts[j+1][i]

            newHorts[j][i] += const * (dv_left + dv_right - 2 * v_net)
            newHorts[j][i+1] += const * (dv_left + dv_right - 2 * v_net)

    return newHorts, newVerts
                
def solve_convection(horts, verts, temps):
    newTemps = temps.copy()

    const = PARAMS["T_STEP"] / PARAMS["L_SCALE"]
    
    for i in range(PARAMS["HEIGHT"]):
        for j in range(PARAMS["WIDTH"]):
            # regions of solidity don't convect (duh)
            if solids[i][j] == 0:
                continue

            # square one above
            if horts[j][i] > 0: # flow is into square (positive side, positive flow)
                # if not along top edge, grab from square above
                if i > 0:
                    top_t = temps[i-1][j]

                # otherwise, check BC for top
                else:
                    if EDGES["TOP"]["TYPE"] == "WALL":
                        top_t = 0

            # if the flow through the top edge is out of the square, use the square's temp
            else:
                top_t = temps[i][j]

            # square one down
            if horts[j][i+1] < 0: # flow is into square (negative side, negative flow)
                if i < PARAMS["HEIGHT"] - 1:
                    bot_t = temps[i+1][j]
                    
                else:
                    if EDGES["BOTTOM"]["TYPE"] == "WALL":
                        top_t = 0
                        
            else:
                bot_t = temps[i][j]


            # square one to left
            if verts[i][j] > 0: # flow is into square (positive side, positive flow)
                if j > 0:
                    lef_t = temps[i][j-1]
                    
                else:
                    if EDGES["LEFT"]["TYPE"] == "WALL":
                        top_t = 0

            else:
                lef_t = temps[i][j]
            
            # square one to right
            if verts[i][j+1] < 0: # flow is into square (negative side, negative flow)
                if j < PARAMS["WIDTH"] - 1:
                    rig_t = temps[i][j+1]
                    
                else:
                    if EDGES["RIGHT"]["TYPE"] == "WALL":
                        top_t = 0

            else:
                rig_t = temps[i][j]

            newTemps[i][j] += (top_t * horts[j][i] + lef_t * verts[i][j] - bot_t * horts[j][i+1] - rig_t * verts[i][j+1]) * const

    return newTemps

def get_vorticities(horts, verts):
    vorts = []
    for i in range(PARAMS["HEIGHT"]):
        vorts.append([])
        for j in range(PARAMS["WIDTH"]):
            vorticity = 0

            # dv/dx
            if j > 0:
                vorticity += (horts[i][j] - horts[i][j-1]) / PARAMS["L_SCALE"]
            if j < WIDTH - 1:
                vorticity += (horts[i][j+1] - horts[i][j]) / PARAMS["L_SCALE"]

            # du/dy
            if i > 0:
                vorticity += -(verts[j][i] - verts[j][i-1]) / PARAMS["L_SCALE"]
            if i < HEIGHT - 1:
                vorticity += -(verts[j][i+1] - verts[j][i1]) / PARAMS["L_SCALE"]
            
            vorts[i].append(vorticity)

    return vorts

def generateImage(temps, minT, maxT, pressures, minP, maxP, resid, t_elapsed, frames):
    newImg = Image.new(mode='RGB', size=(int(PARAMS["WIDTH"]), int(PARAMS["HEIGHT"] * 1.5)), color=WHITE) # temps
    newImg2 = Image.new(mode='RGB', size=(int(PARAMS["WIDTH"]), int(PARAMS["HEIGHT"] * 1.5)), color=WHITE) # pressures

    for i in range(PARAMS["HEIGHT"]):
        for j in range(PARAMS["WIDTH"]):
            try:
                newImg.putpixel((j,i), temp_to_color(temps[i][j], minT, maxT))
            except:
                newImg.putpixel((j,i), (0,0,0))
                keep_going = False

            if solids[i][j] == 0:
                newImg2.putpixel((j,i), SOLID_COLOR)

            else:
                try:
                    newImg2.putpixel((j,i), pres_to_color(pressures[i][j], minP, maxP))
                except:
                    newImg2.putpixel((j,i), (255,0,0))
                    keep_going = False

    d = ImageDraw.Draw(newImg)
    d2 = ImageDraw.Draw(newImg2)

    d.text((PARAMS["WIDTH"] * 0.01, PARAMS["HEIGHT"] * 1.1), "t="+str(t_elapsed), fill=BLACK)
    #d.text((PARAMS["WIDTH"] * 0.4, PARAMS["HEIGHT"] * 1.1), "resid "+str(resid), fill=BLACK)    
    d.text((PARAMS["WIDTH"] * 0.01, PARAMS["HEIGHT"] * 1.3), "max T "+str(maxT), fill=BLACK)
    d.text((PARAMS["WIDTH"] * 0.01, PARAMS["HEIGHT"] * 1.2), "min T "+str(minT), fill=BLACK)

    d2.text((PARAMS["WIDTH"] * 0.01, PARAMS["HEIGHT"] * 1.1), "t="+str(t_elapsed), fill=BLACK)
    #d2.text((PARAMS["WIDTH"] * 0.4, PARAMS["HEIGHT"] * 1.1), "resid "+str(resid), fill=BLACK)
    d2.text((PARAMS["WIDTH"] * 0.01, PARAMS["HEIGHT"] * 1.3), "max P "+str(maxP), fill=BLACK)
    d2.text((PARAMS["WIDTH"] * 0.01, PARAMS["HEIGHT"] * 1.2), "min P "+str(minP), fill=BLACK)

    frameName = str(frames)

    while len(frameName) < len(str(MAX_FRAMES)) - 1:
        frameName = "0" + frameName

    newImg.save("temp_frame" + str(frameName) + ".jpeg", "JPEG")
    newImg2.save("pres_frame" + str(frameName) + ".jpeg", "JPEG")

def generate_vorts_image(horts, verts, resid, t_elapsed):
    newImg = Image.new(mode='RGB', size=(int(PARAMS["WIDTH"]), int(PARAMS["HEIGHT"] * 1.5)), color=WHITE)

    vorticities = get_vorticities(horts, verts)

    minVort, maxVort = getBounds(vorticities)

    for i in range(PARAMS["HEIGHT"]):
        for j in range(PARAMS["WIDTH"]):
            try:
                newImg.putpixel((j,i), vort_to_color(vorticities[i][j], minVort, maxVort))
            except:
                newImg.putpixel((j,i), (128,0,0))

    d = ImageDraw.Draw(newImg)

    d.text((PARAMS["WIDTH"] * 0.01, PARAMS["HEIGHT"] * 1.1), "t="+str(t_elapsed), fill=BLACK)
    #d.text((PARAMS["WIDTH"] * 0.5, PARAMS["HEIGHT"] * 1.1), "resid "+str(resid), fill=BLACK)    
    d.text((PARAMS["WIDTH"] * 0.01, PARAMS["HEIGHT"] * 1.3), "max w "+str(maxVort), fill=BLACK)
    d.text((PARAMS["WIDTH"] * 0.01, PARAMS["HEIGHT"] * 1.2), "min w "+str(minVort), fill=BLACK)

    frameName = str(frames)

    while len(frameName) < len(str(MAX_FRAMES)) - 1:
        frameName = "0" + frameName

    newImg.save("vort_frame"+str(frameName)+".jpeg", "JPEG")

def vort_to_color(vort, minVort, maxVort):
    if abs(minVort) > abs(maxVort):
        k = -minVort
    else:
        k = maxVort

    return int(255 / (1 + math.e ^ (6 * vort / k)))

def getBounds(array):
    minimum = array[0][0]
    maximum = array[0][0]

    for i in range(len(array)):
        for j in range(len(array[0])):
            if array[i][j] < minimum:
                minimum = array[i][j]

            elif array[i][j] > maximum:
                maximum = array[i][j]

    return minimum, maximum

def init_params():
    fileName = input("File name (include extension): ")

    try:
        inFile = open(fileName, 'r')
    except:
        print("No such file")

    params = {} # system parameters
    models = {} # empty dict
    edges = {} # boundary conditions
    properties = {} # fluid (and possibly solid) properties
    fluxes = [] # heat fluxes
    ports = [] # inlets/outlets that don't cover an entire edge
    conds = [] # initial conditions within the domain

    inFileLines = inFile.readlines()

    for line in inFileLines:
        splitLine = line.rstrip("\n").split(" ")


        header = splitLine[0].upper()

        if header == "WIDTH":
            params[header] = int(splitLine[1])

        elif header == "HEIGHT":
            params[header] = int(splitLine[1])

        elif header == "L_SCALE":
            params[header] = float(splitLine[1])

        elif header == "T_STEP":
            try:
                params[header] = float(splitLine[1])
            except:
                params[header] = frac_to_float(splitLine[1])

        elif header == "ERR_THRESH":
            params[header] = float(splitLine[1])

        elif header == "MAX_ITERS":
            params[header] = int(splitLine[1])

        elif header == "G":
            params[header] = float(splitLine[1])

        elif header == "ANIM_FREQ":
            params[header] = int(splitLine[1])

        elif header == "SAVE_FREQ":
            params[header] = int(splitLine[1])

        elif header == "THERM_DIFFUSE":
            models[header] = bool(int(splitLine[1]))
            params["TEMP_ITERS"] = int(splitLine[2])

        elif header == "MOMENTUM":
            models[header] = bool(int(splitLine[1]))

            if models[header]:
                params["RELAX"] = float(splitLine[3])

        elif header == "CONVECTION":
            models[header] = bool(int(splitLine[1]))

        elif header == "BUOYANCY":
            models[header] = bool(int(splitLine[1]))

        elif header == "VISCOUS":
            models[header] = bool(int(splitLine[1]))

        elif header in ("FLUID","SOLID"):
            properties[header] = {}
            for token in range(1,len(splitLine),2):
                properties[header][splitLine[token]] = float(splitLine[token+1])

        elif header in ("TOP", "BOTTOM", "LEFT", "RIGHT"):
            edges[header] = {}

            edgeType = splitLine[1]

            edges[header]["TYPE"] = edgeType

            if edgeType == "WALL":
                if splitLine[2] == "AD":
                    edges[header]["AD"] = True
                    
                elif splitLine[2] == "FLUX":
                    edges[header]["AD"] = False
                    edges[header]["FLUX"] = True
                    edges[header]["TEMP"] = False
                    edges[header]["QUANT"] = float(splitLine[3])

                elif splitLine[2] == "TEMP":
                    edges[header]["AD"] = False
                    edges[header]["FLUX"] = False
                    edges[header]["TEMP"] = True
                    edges[header]["QUANT"] = float(splitLine[3])

        elif header == "FLUX":
            newFlux = {}

            newFlux["EDGE"] = splitLine[1]

            newFlux["START"] = int(splitLine[2])
            newFlux["END"] = int(splitLine[3])

            newFlux["QUANT"] = float(splitLine[4])

            fluxes.append(newFlux)

        elif header == "COND":
            newCond = {"TYPE": splitLine[1]}
            newCond["QUANT"] = float(splitLine[2])
            newCond["SHAPE"] = splitLine[3]
            newCond["MEASURE"] = splitLine[4:]

            for cond in newCond["MEASURE"]:
                cond = int(cond)

            conds.append(newCond)

    return params, models, edges, properties, fluxes, ports, conds

def frac_to_float(frac):
    numer, denom = frac.split("/")

    numer = int(numer)
    denom = int(denom)

    return numer / denom

# function to save data arrays into a readable format
def save_data(fileName, labels, dims, arrays):
    writeFile = open(fileName, mode='wt', encoding='ascii')

    output = ""

    for x in range(len(labels)):
        output = output + labels[x]+" "+str(dims[x])+"\n"

        if dims[x] == 1:
            for i in range(len(arrays[x])):
                output = output + join_nums(arrays[x], ",") +"\n\n"

        elif dims[x] == 2:
            for i in range(len(arrays[x])):
                output = output + join_nums(arrays[x][i],",") +"\n"
            output = output + "\n"


    #print(output)

    writeFile.write(output)

    writeFile.close()

def read_data(fileName):
    inFile = open(fileName, 'r')

    lines = inFile.readlines()

    curData = ""
    dataHold = {}

    for line in lines:
        if curData == "":
            if line == "" or line == "\n":
                continue

            else:
                splitLine = line.rstrip("\n").split(" ")

                if len(splitLine) > 1:
                    curData = splitLine[0]
                    print(curData)
                    numDims = int(splitLine[1])

                    if numDims > 1:
                        dataHold[curData] = []
                        lineNo = 0

        else:
            if numDims == 1:
                dataHold[curData] = line.rstrip("\n").split(",")

                for a in range(len(dataHold[curData])):
                    dataHold[curData][a] = float(dataHold[curData][a])

                curData = ""

            elif numDims == 2:
                if line == "\n" or line =="":
                    curData = ""

                else:
                
                    dataHold[curData].append(line.split(","))

                    for a in range(len(dataHold[curData][lineNo])):
                        try:
                            dataHold[curData][lineNo][a] = float(dataHold[curData][lineNo][a])
                        except:
                            b = 1

                    lineNo += 1

    inFile.close()

    return dataHold

def join_nums(array, sym):
    output = ""

    for x in range(len(array)):
        output = output + str(array[x])

        if x < len(array) - 1:
            output = output + sym

    return output

# call main function
main()
