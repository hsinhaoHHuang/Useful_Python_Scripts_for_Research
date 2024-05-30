import numpy as np

# =============================================================================================================
# Table of Content
# =============================================================================================================
# Step 01. Set precision to output
# Step 02. Read the text file of particle initial condition
# Step 03. Set the arrays for output PAR_IC
# Step 04. Print information for output PAR_IC
# Step 05. Write the output PAR_IC to the file
# =============================================================================================================


# =============================================================================================================
# Step 01. Set precision to output
Float8           = False       # whether the PAR_IC is in double precision

if Float8:
    PAR_IC_dtype = np.double
else:
    PAR_IC_dtype = np.single

# =============================================================================================================

# =============================================================================================================
# Step 02. Read the text file of particle initial condition
Particle_text_filename  = "Particle_000000.txt" # a file

print("")
print("Loading %s ... "%Particle_text_filename)

with open(Particle_text_filename, 'r') as f:
    for line in f:
        if line.startswith('#'):
            header = line
        elif line.startswith( '\n' ):
            continue
        else:
            break #stop when there are no more #

    f.close()

Particle_text_header = header[1:].strip().split()
Particle_text_table  = np.genfromtxt( Particle_text_filename, delimiter=None, comments='#',
                                      names=Particle_text_header, dtype=np.double, encoding=None )
print("done!")

# =============================================================================================================

# =============================================================================================================
# Step 03. Set the arrays for output PAR_IC
Array_Input_ParMass  = Particle_text_table["ParMass"]
Array_Input_ParPosX  = Particle_text_table["ParPosX"]
Array_Input_ParPosY  = Particle_text_table["ParPosY"]
Array_Input_ParPosZ  = Particle_text_table["ParPosZ"]
Array_Input_ParVelX  = Particle_text_table["ParVelX"]
Array_Input_ParVelY  = Particle_text_table["ParVelY"]
Array_Input_ParVelZ  = Particle_text_table["ParVelZ"]
Array_Input_ParType  = Particle_text_table["ParType"]

Num_Par              = Array_Input_ParMass.size  # number of particles

Array_Output_ParMass = Array_Input_ParMass.astype( dtype=PAR_IC_dtype )
Array_Output_ParPosX = Array_Input_ParPosX.astype( dtype=PAR_IC_dtype )
Array_Output_ParPosY = Array_Input_ParPosY.astype( dtype=PAR_IC_dtype )
Array_Output_ParPosZ = Array_Input_ParPosZ.astype( dtype=PAR_IC_dtype )
Array_Output_ParVelX = Array_Input_ParVelX.astype( dtype=PAR_IC_dtype )
Array_Output_ParVelY = Array_Input_ParVelY.astype( dtype=PAR_IC_dtype )
Array_Output_ParVelZ = Array_Input_ParVelZ.astype( dtype=PAR_IC_dtype )
Array_Output_ParType = Array_Input_ParType.astype( dtype=PAR_IC_dtype )

# =============================================================================================================

# =============================================================================================================
# Step 04. Print information for output PAR_IC
print("")
print(f"{Num_Par          = }")
print("")
print(f"{Float8           = }")
print("")
print("Information of the first particle:")
print("")
print(f'{Array_Input_ParMass[0]  =: 24.16e}')
print(f'{Array_Input_ParPosX[0]  =: 24.16e}')
print(f'{Array_Input_ParPosY[0]  =: 24.16e}')
print(f'{Array_Input_ParPosZ[0]  =: 24.16e}')
print(f'{Array_Input_ParVelX[0]  =: 24.16e}')
print(f'{Array_Input_ParVelY[0]  =: 24.16e}')
print(f'{Array_Input_ParVelZ[0]  =: 24.16e}')
print(f'{Array_Input_ParType[0]  =: 24.16e}')
print("")

# =============================================================================================================

# =============================================================================================================
# Step 05. Write the output PAR_IC to the file
output_file = "PAR_IC"

print("")
print("Writing PAR_IC to output file...")
with open(output_file, 'wb') as f:

    f.write( Array_Output_ParMass.tobytes() )
    f.write( Array_Output_ParPosX.tobytes() )
    f.write( Array_Output_ParPosY.tobytes() )
    f.write( Array_Output_ParPosZ.tobytes() )
    f.write( Array_Output_ParVelX.tobytes() )
    f.write( Array_Output_ParVelY.tobytes() )
    f.write( Array_Output_ParVelZ.tobytes() )
    f.write( Array_Output_ParType.tobytes() )

    f.close()

print("done!")

# =============================================================================================================
