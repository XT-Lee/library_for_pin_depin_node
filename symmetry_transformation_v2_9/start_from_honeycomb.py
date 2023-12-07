import hoomd
import hoomd.md
from hoomd import azplugins
import numpy
#import ex_render
import hoomd.azplugins.multi_positions as mp
from matplotlib import pyplot
import hoomd.azplugins.sequence_generator as sg

a_set=3*0.816
seq=sg.sequence()
seq.generate_honeycomb(a=a_set,n=[8,12])

snap=seq.system.take_snapshot()
box=snap.box
points=snap.particles.position[:]
