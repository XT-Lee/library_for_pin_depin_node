
import numpy
#interaction
epsilon=300
kappa=0.25
#substrate
k=800.0
rcut=15.0

energy_interaction=epsilon/rcut*numpy.exp(-kappa*rcut)

#pt=numpy.exp(1)
print('energy_thermal     1')
print('energy_interaction '+str(energy_interaction.astype(int)))
print('energy_interaction '+str(energy_interaction))
