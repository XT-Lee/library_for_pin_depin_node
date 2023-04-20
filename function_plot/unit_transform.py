class unit_transform_module:
    def __init__(self,energy,length,mass):
        R"""
            import function_plot.unit_transform as ut
            stm = ut.unit_transform_module(1,1,1)
            stm.gen_si_unit()
            print(stm.length)
        """
        self.energy = energy
        self.length = length
        self.mass = mass

    def gen_si_unit(self):
        self.si.energy = 1
        self.si.energy_name = 'J'

class si_unit:
    def __init__(self):
        self.energy = 1
        self.length = 1
        self.energy_name = 'J'
        self.length_name = 'm'
        self.mass_name = 'kg'
        