import os
import argparse
import numpy as np

# =============================================================================================================
# Table of Content
# =============================================================================================================
# Step 01. Set parameters in Input__TestProb
# Step 02. Set parameters in Input__Parameter
# Step 03. Load the Input__UM_IC_RefineRegion
# Step 04. Set parameters for output UM_IC
# Step 05. Print the UM_IC information
# Step 06. Define the functions for data construction
# Step 07. Construct the output UM_IC
# Step 08. Write the output UM_IC to the file
# =============================================================================================================


# =============================================================================================================
# Step 01. Set parameters in Input__TestProb
Plummer_Rho0     =  1.0                 # peak density
Plummer_R0       =  0.1                 # scale radius
Plummer_Center   = [1.50, 1.50, 1.50]   # central coordinates
Plummer_BulkVel  = [0.00, 0.00, 0.00]   # bulk velocity
Plummer_GasMFrac =  0.25                # gas mass fraction

# =======================================================================================

# =============================================================================================================
# Step 02. Set parameters in Input__Parameter
GAMMA            = 1.666666667 # ratio of specific heats (i.e., adiabatic index)
NEWTON_G         = 1.0         # gravitational constant

UM_IC_BoxSize_x  = 3.0         # box size of the UM_IC in the x direction
UM_IC_N_x_base   = 64          # number of base-level cells in the UM_IC in the x direction
UM_IC_N_y_base   = 64          # ...                                            y direction
UM_IC_N_z_base   = 64          # ...                                            z direction

MAX_LEVEL        = 2           # maximum refinement level

PatchSize        = 8           # PATCH_SIZE in GAMER
Float8           = False       # whether the UM_IC is in double precision

if Float8:
    UM_IC_dtype = np.double
else:
    UM_IC_dtype = np.single

# =======================================================================================

# =============================================================================================================
# Step 03. Load the Input__UM_IC_RefineRegion
Input__UM_IC_RefineRegion_filename  = "Input__UM_IC_RefineRegion" # a file to specify the refinement region of the input multi-level UM_IC

print("")
print("Loading %s ... "%Input__UM_IC_RefineRegion_filename)

with open(Input__UM_IC_RefineRegion_filename, 'r') as f:
    for line in f:
        if line.startswith('#'):
            header = line
        else:
            break #stop when there are no more #

    f.close()

Input__UM_IC_RefineRegion_header = header[1:].strip().split()
Input__UM_IC_RefineRegion_table  = np.genfromtxt( Input__UM_IC_RefineRegion_filename, delimiter=None, comments='#',
                                                  names=Input__UM_IC_RefineRegion_header, dtype=None, encoding=None )

print("done!")

# =============================================================================================================

# =============================================================================================================
# Step 04. Set parameters for output UM_IC

UM_IC_dh_base     = UM_IC_BoxSize_x/UM_IC_N_x_base # base-level cell size of the input UM_IC
UM_IC_BoxSize_y   = UM_IC_dh_base*UM_IC_N_y_base   # box size of the input UM_IC in the y direction
UM_IC_BoxSize_z   = UM_IC_dh_base*UM_IC_N_z_base   # ...                                z direction


UM_IC_NP_Skip_xL  = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # number of patches on the parent level to be skipped in the x direction from the left edge of the parent refinement region
UM_IC_NP_Skip_xR  = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # ...                                                        x ...                right ...
UM_IC_NP_Skip_yL  = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # ...                                                        y ...                left  ...
UM_IC_NP_Skip_yR  = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # ...                                                        y ...                right ...
UM_IC_NP_Skip_zL  = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # ...                                                        z ...                left ...
UM_IC_NP_Skip_zR  = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # ...                                                        z ...                right ...
UM_IC_NP_x        = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # number of patches on each level in the x direction
UM_IC_NP_y        = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # ...                                    y direction
UM_IC_NP_z        = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # ...                                    z direction
UM_IC_N_x         = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # number of cells on each level in the x direction
UM_IC_N_y         = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # ...                                  y direction
UM_IC_N_z         = np.zeros( MAX_LEVEL+1, dtype=np.uint32 )  # ...                                  z direction
UM_IC_x0          = np.zeros( MAX_LEVEL+1 )                   # left edge of the refinement region for each level in the x direction
UM_IC_y0          = np.zeros( MAX_LEVEL+1 )                   # left ...                                                 y direction
UM_IC_z0          = np.zeros( MAX_LEVEL+1 )                   # left ...                                                 z direction
UM_IC_x1          = np.zeros( MAX_LEVEL+1 )                   # right ...                                                x direction
UM_IC_y1          = np.zeros( MAX_LEVEL+1 )                   # right ...                                                y direction
UM_IC_z1          = np.zeros( MAX_LEVEL+1 )                   # right ...                                                z direction
UM_IC_dh          = np.zeros( MAX_LEVEL+1 )                   # cell size of the input UM_IC for each level

# Loop for each level to set the UM_IC structure information
for lv in range(0, MAX_LEVEL+1, 1):

    UM_IC_NP_Skip_xL[lv] = 0                        if lv == 0 else Input__UM_IC_RefineRegion_table["NP_Skip_xL"][(Input__UM_IC_RefineRegion_table["dLv"] == lv)]
    UM_IC_NP_Skip_xR[lv] = 0                        if lv == 0 else Input__UM_IC_RefineRegion_table["NP_Skip_xR"][(Input__UM_IC_RefineRegion_table["dLv"] == lv)]
    UM_IC_NP_Skip_yL[lv] = 0                        if lv == 0 else Input__UM_IC_RefineRegion_table["NP_Skip_yL"][(Input__UM_IC_RefineRegion_table["dLv"] == lv)]
    UM_IC_NP_Skip_yR[lv] = 0                        if lv == 0 else Input__UM_IC_RefineRegion_table["NP_Skip_yR"][(Input__UM_IC_RefineRegion_table["dLv"] == lv)]
    UM_IC_NP_Skip_zL[lv] = 0                        if lv == 0 else Input__UM_IC_RefineRegion_table["NP_Skip_zL"][(Input__UM_IC_RefineRegion_table["dLv"] == lv)]
    UM_IC_NP_Skip_zR[lv] = 0                        if lv == 0 else Input__UM_IC_RefineRegion_table["NP_Skip_zR"][(Input__UM_IC_RefineRegion_table["dLv"] == lv)]
    UM_IC_NP_x      [lv] = UM_IC_N_x_base/PatchSize if lv == 0 else 2*(UM_IC_NP_x[lv-1]-UM_IC_NP_Skip_xL[lv]-UM_IC_NP_Skip_xR[lv])
    UM_IC_NP_y      [lv] = UM_IC_N_y_base/PatchSize if lv == 0 else 2*(UM_IC_NP_y[lv-1]-UM_IC_NP_Skip_yL[lv]-UM_IC_NP_Skip_yR[lv])
    UM_IC_NP_z      [lv] = UM_IC_N_z_base/PatchSize if lv == 0 else 2*(UM_IC_NP_z[lv-1]-UM_IC_NP_Skip_zL[lv]-UM_IC_NP_Skip_zR[lv])
    UM_IC_N_x       [lv] = UM_IC_N_x_base           if lv == 0 else UM_IC_NP_x[lv]*PatchSize
    UM_IC_N_y       [lv] = UM_IC_N_y_base           if lv == 0 else UM_IC_NP_y[lv]*PatchSize
    UM_IC_N_z       [lv] = UM_IC_N_z_base           if lv == 0 else UM_IC_NP_z[lv]*PatchSize
    UM_IC_x0        [lv] = 0.0                      if lv == 0 else UM_IC_x0[lv-1]+UM_IC_NP_Skip_xL[lv]*PatchSize*UM_IC_dh[lv-1]
    UM_IC_y0        [lv] = 0.0                      if lv == 0 else UM_IC_y0[lv-1]+UM_IC_NP_Skip_yL[lv]*PatchSize*UM_IC_dh[lv-1]
    UM_IC_z0        [lv] = 0.0                      if lv == 0 else UM_IC_z0[lv-1]+UM_IC_NP_Skip_zL[lv]*PatchSize*UM_IC_dh[lv-1]
    UM_IC_x1        [lv] = UM_IC_BoxSize_x          if lv == 0 else UM_IC_x1[lv-1]-UM_IC_NP_Skip_xR[lv]*PatchSize*UM_IC_dh[lv-1]
    UM_IC_y1        [lv] = UM_IC_BoxSize_y          if lv == 0 else UM_IC_y1[lv-1]-UM_IC_NP_Skip_yR[lv]*PatchSize*UM_IC_dh[lv-1]
    UM_IC_z1        [lv] = UM_IC_BoxSize_z          if lv == 0 else UM_IC_z1[lv-1]-UM_IC_NP_Skip_zR[lv]*PatchSize*UM_IC_dh[lv-1]
    UM_IC_dh        [lv] = UM_IC_dh_base            if lv == 0 else UM_IC_dh_base/(2**lv)

# =============================================================================================================

# =============================================================================================================
# Step 05. Print the UM_IC information
print("")
print("------------------------------------------------------------------------------------------------")
print("UM_IC information")
print("------------------------------------------------------------------------------------------------")
print(f"{PatchSize         = }")
print(f"{Float8            = }")
print(f"{MAX_LEVEL         = }")
print("")
print(f"{UM_IC_BoxSize_x   = }")
print(f"{UM_IC_BoxSize_y   = }")
print(f"{UM_IC_BoxSize_z   = }")
print("")
print(f"{UM_IC_N_x_base    = }")
print(f"{UM_IC_N_y_base    = }")
print(f"{UM_IC_N_z_base    = }")
print("")
print(f"{UM_IC_dh_base     = }")
print("")
print(f"{UM_IC_NP_Skip_xL  = }")
print(f"{UM_IC_NP_Skip_xR  = }")
print(f"{UM_IC_x0          = }")
print(f"{UM_IC_x1          = }")
print("")
print(f"{UM_IC_NP_Skip_yL  = }")
print(f"{UM_IC_NP_Skip_yR  = }")
print(f"{UM_IC_y0          = }")
print(f"{UM_IC_y1          = }")
print("")
print(f"{UM_IC_NP_Skip_zL  = }")
print(f"{UM_IC_NP_Skip_zR  = }")
print(f"{UM_IC_z0          = }")
print(f"{UM_IC_z1          = }")
print("")
print(f"{UM_IC_NP_x        = }")
print(f"{UM_IC_NP_y        = }")
print(f"{UM_IC_NP_z        = }")
print("")
print(f"{UM_IC_N_x         = }")
print(f"{UM_IC_N_y         = }")
print(f"{UM_IC_N_z         = }")
print("")
print(f"{UM_IC_dh          = }")
print("------------------------------------------------------------------------------------------------")
print("")

# =============================================================================================================

# =======================================================================================
# Step 06. Define the functions for data construction
# Eint from Dens and Pres
def EoS_DensPres2Eint_Gamma( Dens, Pres ):

    Eint = Pres * 1.0 / ( GAMMA - 1.0 )

    return Eint

# Etot from Con and Eint
def Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, Emag ):

    Etot  = 0.5*( (MomX)**2 + (MomY)**2 + (MomZ)**2 ) / Dens
    Etot += Eint
    Etot += Emag

    return Etot

# Set the fluid data
def SetGridIC( x, y, z, Time ):

    # gas share the same density profile as particles (except for different total masses)
    TotM       = 4.0/3.0*np.pi*(Plummer_R0**3)*Plummer_Rho0
    GasRho0    = Plummer_Rho0*Plummer_GasMFrac
    PresBg     = 0.0                                       # background pressure (set to 0.0 by default)

    r2         = (x-Plummer_Center[0])**2 + (y-Plummer_Center[1])**2 + (z-Plummer_Center[2])**2
    a2         = r2 / (Plummer_R0)**2
    Dens1Cloud = GasRho0 * ( 1.0 + a2 )**(-2.5)

    # set fluid variables
    Dens = Dens1Cloud
    MomX = Dens1Cloud*Plummer_BulkVel[0]
    MomY = Dens1Cloud*Plummer_BulkVel[1]
    MomZ = Dens1Cloud*Plummer_BulkVel[2]
    Pres = NEWTON_G*TotM*GasRho0 / ( 6.0*Plummer_R0*(1.0 + a2)**3 ) + PresBg

    # compute the total gas energy
    Eint = EoS_DensPres2Eint_Gamma( Dens, Pres )                    # assuming EoS requires no passive scalars
    Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 )  # do NOT include magnetic energy here

    # set the output
    DENS = Dens
    MOMX = MomX
    MOMY = MomY
    MOMZ = MomZ
    ENGY = Etot

    return DENS, MOMX, MOMY, MOMZ, ENGY

# =============================================================================================================

# =============================================================================================================
# Step 07. Construct the output UM_IC
print("")
print("Constructing the output UM_IC ...")

for lv in range(0, MAX_LEVEL+1, 1):
    print("    lv %d ..."%lv)

    Nx = UM_IC_N_x[lv]
    Ny = UM_IC_N_y[lv]
    Nz = UM_IC_N_z[lv]
    x0 = UM_IC_x0[lv]
    y0 = UM_IC_y0[lv]
    z0 = UM_IC_z0[lv]
    dh = UM_IC_dh[lv]

    UM_IC_thisLv = np.zeros( (5, Nz, Ny, Nx), dtype=UM_IC_dtype )
    Array_Dens   = np.zeros(    (Nz, Ny, Nx), dtype=UM_IC_dtype )
    Array_Momx   = np.zeros(    (Nz, Ny, Nx), dtype=UM_IC_dtype )
    Array_Momy   = np.zeros(    (Nz, Ny, Nx), dtype=UM_IC_dtype )
    Array_Momz   = np.zeros(    (Nz, Ny, Nx), dtype=UM_IC_dtype )
    Array_Engy   = np.zeros(    (Nz, Ny, Nx), dtype=UM_IC_dtype )

    for i in range(0, Nx, 1):
        for j in range(0, Ny, 1):
            for k in range(0, Nz, 1):

                x = x0 + (i+0.5)*dh
                y = y0 + (j+0.5)*dh
                z = z0 + (k+0.5)*dh

                Dens, Momx, Momy, Momz, Engy = SetGridIC( x, y, z, 0.0 )

                Array_Dens[k, j, i] = Dens
                Array_Momx[k, j, i] = Momx
                Array_Momy[k, j, i] = Momy
                Array_Momz[k, j, i] = Momz
                Array_Engy[k, j, i] = Engy

    UM_IC_thisLv[0, :, :, :] = Array_Dens
    UM_IC_thisLv[1, :, :, :] = Array_Momx
    UM_IC_thisLv[2, :, :, :] = Array_Momy
    UM_IC_thisLv[3, :, :, :] = Array_Momz
    UM_IC_thisLv[4, :, :, :] = Array_Engy

    # =============================================================================================================

    # =============================================================================================================
    # Step 08. Write the output UM_IC to the file
    output_file    = "UM_IC"
    output_file_Lv = output_file+"_Lv%02d"%lv

    print("    Writing UM_IC to output file...")
    with open(output_file_Lv, "wb") as f:
        UM_IC_thisLv.tofile(f)
        f.close()

    os.system("cat %s >> %s"%(output_file_Lv, output_file))
    print("    done!")

    # =============================================================================================================

print("done!")

# =============================================================================================================
