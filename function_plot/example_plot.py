import numpy as np
import matplotlib.pyplot as plt
class functions_plot_module:
    R"""
    exp:
        import function_plot.example_plot as ex
        ef = ex.functions_plot_module()
        x,y = ef.generate_power_law()
        ef.plot_function_semilogy(x,y)
        x,y = ef.generate_cumulate(x,y)
        ef.plot_function(x,y)
    """
    def __init__(self) :
        pass

    def generate_power_law(self):
        R"""
        introduction: 
            ln(y) = k*ln(x), linear in loglog,
            when -k > 1 , the cumulation is convergent series
        return: 
            x,y
        """
        x = np.linspace(0.1,10,100)
        y = np.power(x,-2)
        return x,y
    
    def generate_cumulate(self,x,y):
        R"""
        return x,y
        """
        y = tuple(y)
        sumy = np.array(y)
        sz = np.shape(sumy)
        for i in range(sz[0]):
            if i>0:
                sumy[i] = sumy[i]+sumy[i-1]
        return x,sumy

    def generate_exp(self):
        R"""
        introduction: 
            ln(y) = x*ln(e), linear in semilogy
            when e > 1 , the cumulation is convergent series
        return: 
            x,y
        """
        x = np.linspace(0.1,0.9,10)
        y = np.exp(-x)*2/1.1
        return x,y
    
    def generate_gaussian(self,epsilon = 300,sigma = 1,rcut = 2):
        R"""
        introduction: 
            
        return: 
            x,y
        example:
            import function_plot.example_plot as fp
            fpm = fp.functions_plot_module()
            fpm.generate_gaussian(1,1,1)
            fpm.generate_gaussian(10,1,1)
            fpm.generate_gaussian(1,0.6,2)
            fpm.generate_gaussian(10,0.6,2)
        """
        # = 
        
        x = np.linspace(-rcut,rcut,41)
        y = -np.exp(-0.5*((x/sigma)**2))*epsilon + np.exp(-0.5*((rcut/sigma)**2))*epsilon
        delta_energy_per_epsilon = (1 - np.exp(-0.5*((rcut/sigma)**2)))*epsilon
        print('sigma',sigma,'rcut',rcut,'dE',delta_energy_per_epsilon)
        return x,y
    
    def generate_harmonic(self,epsilon = 300):
        R"""
        introduction: 
            ln(y) = x*ln(e), linear in semilogy
            when e > 1 , the cumulation is convergent series
        return: 
            x,y
        """
        # = 
        #sigma = 1
        x = np.linspace(-1,1,21)
        y = 0.5*(x**2)*epsilon - 0.5*epsilon
        return x,y

    def generate_gaussian_force(self,epsilon = 300,sigma = 1):
        R"""
        introduction: 
            ln(y) = x*ln(e), linear in semilogy
            when e > 1 , the cumulation is convergent series
        return: 
            x,y
        """
        x = np.linspace(-1,1,21)
        y = x*np.exp(-0.5*((x/sigma)**2))*epsilon/(sigma**2)
        return x,y
    
    def generate_harmonic_force(self,epsilon = 300):
        R"""
        introduction: 
            ln(y) = x*ln(e), linear in semilogy
            when e > 1 , the cumulation is convergent series
        return: 
            x,y
        """
        # = 
        #sigma = 1
        x = np.linspace(-1,1,21)
        y = x*epsilon
        return x,y

    def generate_kauzmann_glass(self,b,m):
        R"""
        introduction: 
            ln(y) = x*ln(e), linear in semilogy
            when e > 1 , the cumulation is convergent series
        return: 
            x,y
        """
        xm = np.linspace(0.01,1,99)#Tg/T
        x = 1/xm#T/Tg
        y = b/(x-m)#m=T0/Tg
        return x,y
    
    def generate_glass_list(self):
        R"""
        introduction: 
            ln(y) = x*ln(e), linear in semilogy
            when e > 1 , the cumulation is convergent series
        return: 
            x,y
        """
        b=1
        xm = np.linspace(0.01,1,99)#Tg/T
        x = 1/xm#T/Tg
        fig,ax = plt.subplots()
        ax.set_xlabel('$T_g/T$')
        ax.set_ylabel('$\eta$')
        
        for m in np.linspace(0.1,0.9,9):
            y = b/(x-m)#m=T0/Tg
            ax.plot(1/x,y)
        plt.show()

    def generate_compare(self):
        xm = np.linspace(0.01,1,99)#Tg/T
        x = 1/xm#T/Tg
        print(x)

    def plot_function_glass(self,x,y):
        fig,ax = plt.subplots()
        ax.plot(1/x,y)
        ax.set_xlabel('$T_g/T$')
        ax.set_ylabel('$\eta$')
        plt.show()

    def plot_function(self,x,y):
        fig,ax = plt.subplots()
        ax.plot(x,y)
        plt.show()
    
    def plot_function2(self,x,y1,y2):
        R"""
        intro: U/kt = eps ~ k/2; force_gaus < force_harmo,
            when sigma=0.6, the two potential look similar.
        import function_plot.example_plot as ep
        fpm = ep.functions_plot_module()
        eps=300
        x1,y1 = fpm.generate_gaussian(eps,0.6)
        x2,y2 = fpm.generate_harmonic(2*eps)
        fpm.plot_function2(x1,y1,y2)
        """
        fig,ax = plt.subplots()
        ax.plot(x,y1,label='gaus')
        ax.plot(x,y2,label='harm')
        plt.legend()
        plt.show()
    
    def plot_function22(self,x1,y1,x2,y2):
        R"""
        intro: U/kt = eps ~ k/2; force_gaus < force_harmo,
            when epsilon=2k, sigma=0.6,rcut=2, 
            the two potential look similar.
        import function_plot.example_plot as ep
        fpm = ep.functions_plot_module()

        eps=300
        for sigma in [0.6,1]:
            x1,y1 = fpm.generate_gaussian(eps,sigma)
            x2,y2 = fpm.generate_harmonic(2*eps)
            fpm.plot_function22(x1,y1,x2,y2)
            x1,y1 = fpm.generate_gaussian_force(eps,sigma)
            x2,y2 = fpm.generate_harmonic_force(2*eps)
            fpm.plot_function22(x1,y1,x2,y2)
        """
        fig,ax = plt.subplots()
        ax.plot(x1,y1,label='gaus')
        ax.plot(x2,y2,label='harm')
        plt.legend()
        plt.show()
    
    def plot_function_loglog(self,x,y):
        fig,ax = plt.subplots()
        ax.loglog(x,y)
        plt.show()
    
    def plot_function_semilogy(self,x,y):
        fig,ax = plt.subplots()
        ax.semilogy(x,y)
        plt.show()