#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

####################################################
# Genereal Setting
####################################################
dpi        = 150
FONT_SIZE  = 12.0
LINE_WIDTH = 2.0
MARKER_SIZE= 6.0

plt.rcParams['font.size']         = FONT_SIZE
plt.rcParams['figure.titlesize']  = FONT_SIZE
plt.rcParams['figure.dpi']        = dpi

plt.rcParams['axes.titlesize']    = FONT_SIZE
plt.rcParams['axes.labelsize']    = FONT_SIZE
plt.rcParams['axes.labelpad']     = FONT_SIZE/8.0*0.05
plt.rcParams['axes.linewidth']    = LINE_WIDTH

plt.rcParams['legend.fontsize']   = FONT_SIZE
plt.rcParams['lines.linewidth']   = LINE_WIDTH

plt.rcParams['xtick.major.size']  = 4.0*LINE_WIDTH
plt.rcParams['xtick.major.width'] = LINE_WIDTH
plt.rcParams['xtick.minor.size']  = 2.0*LINE_WIDTH
plt.rcParams['xtick.minor.width'] = 0.5*LINE_WIDTH
plt.rcParams['xtick.labelsize']   = FONT_SIZE

plt.rcParams['ytick.major.size']  = 4.0*LINE_WIDTH
plt.rcParams['ytick.major.width'] = LINE_WIDTH
plt.rcParams['ytick.minor.size']  = 2.0*LINE_WIDTH
plt.rcParams['ytick.minor.width'] = 0.5*LINE_WIDTH
plt.rcParams['ytick.labelsize']   = FONT_SIZE

plt.rcParams['font.family']       = 'STIXGeneral'
plt.rcParams['mathtext.fontset']  = 'custom'
plt.rcParams['mathtext.rm']       = 'STIXGeneral:regular'
plt.rcParams['mathtext.it']       = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf']       = 'STIXGeneral:italic:bold'
###################################################################################################

def LoadTable( filename ):
    print( '' )
    print( 'Loading %s ... '%filename )

    with open( filename, 'r' ) as f:
        for line in f:
            if line.startswith( '#' ):
                header = line
            elif line.startswith( '\n' ):
                continue
            else:
                break #stop when there are no more #

        f.close()

    file_header = header[1:].strip().split()
    file_table  = np.genfromtxt( filename, delimiter=None, comments='#',
                                 names=file_header, dtype=None, encoding=None )

    print( 'done!' )
    return file_table

def main() -> None:

    Table1 = LoadTable( 'Record1st/Record__Conservation' )
    Table2 = LoadTable( 'Record2nd/Record__Conservation' )
    Table3 = LoadTable( 'Record3rd/Record__Conservation' )
    Table4 = LoadTable( 'Record4th/Record__Conservation' )

    Mass_Psi_Ref = Table1['Mass_Psi'][0]
    Etot_Psi_Ref = Table1['Etot_Psi'][0]
    StepShow = 1000

    # create figure
    fig = plt.figure( figsize=(8,6) )

    # plot the lines with the colors
    ax = fig.add_subplot(111)

    ax.set_title( 'Conservation Error' )
    ax.plot( Table1['Time'], (Table1['Mass_Psi']-Mass_Psi_Ref)/np.abs(Mass_Psi_Ref)*100.0, '-', color='C0', label='Mass'  )
    ax.plot( Table2['Time'], (Table2['Mass_Psi']-Mass_Psi_Ref)/np.abs(Mass_Psi_Ref)*100.0, '-', color='C0', label='_Mass' )
    ax.plot( Table3['Time'], (Table3['Mass_Psi']-Mass_Psi_Ref)/np.abs(Mass_Psi_Ref)*100.0, '-', color='C0', label='_Mass' )
    ax.plot( Table4['Time'], (Table4['Mass_Psi']-Mass_Psi_Ref)/np.abs(Mass_Psi_Ref)*100.0, '-', color='C0', label='_Mass' )
    ax.plot( Table1['Time'], (Table1['Etot_Psi']-Etot_Psi_Ref)/np.abs(Etot_Psi_Ref)*100.0, '-', color='C1', label='Etot'  )
    ax.plot( Table2['Time'], (Table2['Etot_Psi']-Etot_Psi_Ref)/np.abs(Etot_Psi_Ref)*100.0, '-', color='C1', label='_Etot' )
    ax.plot( Table3['Time'], (Table3['Etot_Psi']-Etot_Psi_Ref)/np.abs(Etot_Psi_Ref)*100.0, '-', color='C1', label='_Etot' )
    ax.plot( Table4['Time'], (Table4['Etot_Psi']-Etot_Psi_Ref)/np.abs(Etot_Psi_Ref)*100.0, '-', color='C1', label='_Etot' )
    ax.axvline( Table1['Time'][(Table1['Step']==StepShow)][0]*code_time_to_Myr, linestyle='--', color='C4', label='Step = %d'%StepShow )
    ax.axhline( Table1['Mass_Psi_RErr'][(Table1['Step']==StepShow)][0]*100.0,   linestyle=':',  color='C0', label='%5.2f%%'%(Table1['Mass_Psi_RErr'][(Table1['Step']==StepShow)][0]*100.0) )
    ax.axhline( Table1['Etot_Psi_RErr'][(Table1['Step']==StepShow)][0]*100.0,   linestyle=':',  color='C1', label='%5.2f%%'%(Table1['Etot_Psi_RErr'][(Table1['Step']==StepShow)][0]*100.0) )
    ax.set_xlabel( 'Time (code_time)' )
    ax.set_ylabel( 'Relative Error (%)' )
    ax.set_xlim( left=0.0 )
    ax.set_ylim( -100.0, 100.0 )
    ax.grid()
    ax.set_yticks([ i for i in range(-100, 100, 10) ])
    ax.xaxis.get_ticklocs( minor=True )
    ax.yaxis.get_ticklocs( minor=True )
    ax.minorticks_on()
    ax.xaxis.set_ticks_position( 'both' )
    ax.yaxis.set_ticks_position( 'both' )
    ax.tick_params( which='both', direction='in' )
    ax.legend()

    # save to file
    plt.tight_layout( pad=0.2 )
    fig.set_dpi(dpi)
    fig.savefig( 'fig_ConservationError.png', dpi=dpi )

    ax.set_xscale('symlog', linthresh=(Table1['Time'][(Table1['Step']==1)][0]*code_time_to_Myr))
    ax.set_yscale('symlog', linthresh=1e-4)
    # save to file
    plt.tight_layout( pad=0.2 )
    fig.set_dpi(dpi)
    fig.savefig( 'fig_ConservationError_log.png', dpi=dpi )

if __name__ == '__main__':
    main()
