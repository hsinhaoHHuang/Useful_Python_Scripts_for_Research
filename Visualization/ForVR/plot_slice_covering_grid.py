import yt
import numpy as np

filename = 'CoveringGrid.h5'
field    = ('grid','Dens')

ds = yt.load( filename )

sz = yt.SlicePlot( ds, 'z', field, center='c' )
sz.save( 'fig_CoveringGrid_Slice_%s.png'%field[1] )
