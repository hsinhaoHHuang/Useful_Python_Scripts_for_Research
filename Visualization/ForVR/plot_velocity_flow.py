import yt
import numpy as np
import add_ELBDM_derived_fields

ds = yt.load( '../Data_000001' )
field         = ('gamer','S')
colormap_S    = 'rainbow'

add_ELBDM_derived_fields.Add_ELBDM_derived_fields( ds )
dh = (ds.domain_width[0]/ds.domain_dimensions[0])/2
N = 1024

sz = yt.SlicePlot( ds, 'z', field, center=[ 320+102.206, 320+23.3333, 320+27.9842 ]*dh, width=10*dh, buff_size=[N, N] )
sz.annotate_quiver( ('gamer','bulk_velocity_x'), ('gamer','bulk_velocity_y'), ('gamer','bulk_velocity_magnitude'), factor=16, cmap='inferno_r' )
sz.annotate_streamlines( ('gamer','bulk_velocity_x'), ('gamer','bulk_velocity_y'))
sz.set_log( field, False )
sz.set_cmap( field, colormap_S )
sz.save()

sy = yt.SlicePlot( ds, 'y', field, center=[ 320+102.206, 320+23.3333, 320+27.9842 ]*dh, width=10*dh, buff_size=[N, N] )
sy.annotate_quiver( ('gamer','bulk_velocity_z'), ('gamer','bulk_velocity_x'), ('gamer','bulk_velocity_magnitude'), factor=16, cmap='inferno_r' )
sy.annotate_streamlines( ('gamer','bulk_velocity_z'), ('gamer','bulk_velocity_x'))
sy.set_log( field, False )
sy.set_cmap( field, colormap_S )
sy.save()

sx = yt.SlicePlot( ds, 'x', field, center=[ 320+102.206, 320+23.3333, 320+27.9842 ]*dh, width=10*dh, buff_size=[N, N] )
sx.annotate_quiver( ('gamer','bulk_velocity_y'), ('gamer','bulk_velocity_z'), ('gamer','bulk_velocity_magnitude'), factor=16, cmap='inferno_r' )
sx.annotate_streamlines( ('gamer','bulk_velocity_y'), ('gamer','bulk_velocity_z'))
sx.set_log( field, False )
sx.set_cmap( field, colormap_S )
sx.save()

sf  = yt.OffAxisSlicePlot( ds,
                           [ 0.670918, -0.232255, 0.70422 ],
                           field,
                           center=[ 320+102.206, 320+23.3333, 320+27.9842 ]*dh,
                           width=10*dh,
                           buff_size=[N, N] )
sf.set_cmap( field, colormap_S )
sf.set_log( field, False )
sf.save()
