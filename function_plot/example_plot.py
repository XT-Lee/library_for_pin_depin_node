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
    
    def plot_function_loglog(self,x,y):
        fig,ax = plt.subplots()
        ax.loglog(x,y)
        plt.show()
    
    def plot_function_semilogy(self,x,y):
        fig,ax = plt.subplots()
        ax.semilogy(x,y)
        plt.show()