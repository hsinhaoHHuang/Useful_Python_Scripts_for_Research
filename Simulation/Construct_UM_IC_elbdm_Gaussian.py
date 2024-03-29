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
Wave   = 1
Hybrid = 2
Scheme = Wave

Gau_v0        = 0.0     #  mean velocity
Gau_Width     = 0.1     #  Gaussian width
Gau_Center    = 0.25    #  Gaussian center
Gau_XYZ       = 0       #  wave propagation direction (0/1/2 --> x/y/z)
Gau_PeriodicN = 10      #  periodic boundary condition
                        #            // (0 = non-periodic, >0 = number of periodic images each side)
ELBDM_ETA     = 6.0

# =======================================================================================

# =============================================================================================================
# Step 02. Set parameters in Input__Parameter

UM_IC_BoxSize_x  = 1.0         # box size of the UM_IC in the x direction
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
# Set the fluid data
def SetGridIC( x, y, z, Time ):

    if Gau_XYZ == 0:
       r = x
       BoxSize = UM_IC_BoxSize_x
    if Gau_XYZ == 1:
       r = y
       BoxSize = UM_IC_BoxSize_y
    if Gau_XYZ == 2:
       r = z
       BoxSize = UM_IC_BoxSize_z

    Gau_Const1 = 1.0 + (  Time / ( ELBDM_ETA*(Gau_Width)**2 ) )**2
    Gau_Theta1 = -0.5*np.arccos( ( Gau_Const1 )**(-0.5)  )
    Real = 0.0
    Imag = 0.0

##  n=0, m=0: original wave packet
##  n>0, m=0/1: images for periodic BC on the plus(+)/minus(-) direction
    for n in range(0, Gau_PeriodicN+1, 1):
        for m in range(0, [2,1][n==0], 1):
            Center     = Gau_Center + n*(1-2*m)*BoxSize
            dr1        = r -     Gau_v0*Time - Center
            dr2        = r - 0.5*Gau_v0*Time - Center
            Gau_Const2 = ( (Gau_Width**2)*np.pi*Gau_Const1 )**(-0.25)*np.exp( -0.5*( dr1/Gau_Width )**2/Gau_Const1 )
            Gau_Theta2 = 0.5*( dr1**2 )*ELBDM_ETA*Time/( ( ELBDM_ETA*(Gau_Width**2) )**2 + (Time)**2  ) + Gau_v0*ELBDM_ETA*dr2

            Real += Gau_Const2*np.cos( Gau_Theta1 + Gau_Theta2 )
            Imag += Gau_Const2*np.sin( Gau_Theta1 + Gau_Theta2 )

    return Real, Imag

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

    UM_IC_thisLv = np.zeros( (2, Nz, Ny, Nx), dtype=UM_IC_dtype )
    Array_Real   = np.zeros(    (Nz, Ny, Nx), dtype=UM_IC_dtype )
    Array_Imag   = np.zeros(    (Nz, Ny, Nx), dtype=UM_IC_dtype )

    for i in range(0, Nx, 1):
        for j in range(0, Ny, 1):
            for k in range(0, Nz, 1):

                x = x0 + (i+0.5)*dh
                y = y0 + (j+0.5)*dh
                z = z0 + (k+0.5)*dh

                Real, Imag  = SetGridIC( x, y, z, 0.0 )
                Array_Real[k, j, i] = Real
                Array_Imag[k, j, i] = Imag

    # Phase and Density
    Array_Phas = np.arctan2( Array_Imag, Array_Real )
    Array_Dens = Array_Real**2 + Array_Imag**2

    # Unwrap phase
    Array_Phas = np.unwrap( Array_Phas, axis=0 )
    Array_Phas = np.unwrap( Array_Phas, axis=1 )
    Array_Phas = np.unwrap( Array_Phas, axis=2 )

    if Scheme == Hybrid:
        UM_IC_thisLv[0, :, :, :] = Array_Dens
        UM_IC_thisLv[1, :, :, :] = Array_Phas

    elif Scheme == Wave:
        UM_IC_thisLv[0, :, :, :] = Array_Real
        UM_IC_thisLv[1, :, :, :] = Array_Imag

    # =============================================================================================================

    # =============================================================================================================
    # Step 08. Write the output UM_IC to the file
    if Scheme == Hybrid:
        output_file = "UM_IC_hybrid"
    elif Scheme == Wave:
        output_file = "UM_IC_wave"

    output_file_Lv = output_file+"_Lv%02d"%lv

    print("        Writing UM_IC to output file...")
    with open(output_file_Lv, "wb") as f:
        UM_IC_thisLv.tofile(f)
        f.close()

    os.system("cat %s >> %s"%(output_file_Lv, output_file))
    print("        done!")
    # =============================================================================================================

    print("    done!")

    # =============================================================================================================

print("done!")

# =============================================================================================================
