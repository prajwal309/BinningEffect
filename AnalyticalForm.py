import sympy
import numpy as np

Log2 = np.log(2)
Log2 = np.log(2)


#Defining all the parameters
gamma_D = sympy.symbols("gamma_{D}", positive=True)
gamma_L = sympy.symbols("gamma_{L}", positive=True)
nu = sympy.symbols("nu", positive=True)
nu_0 = sympy.symbols("nu_0", constant=True)

#Based on the derivation in Liu et. al 2001
#d = sympy.symbols("d", positive=True)
d = (gamma_L - gamma_D)/(gamma_L+gamma_D)

alpha = 0.18121
beta = 0.023665*sympy.exp(0.6*d)+0.00418*sympy.exp(-1.9*d)

gamma_V = (gamma_L + gamma_D)*(1- alpha*(1-d*d) - beta*sympy.sin(sympy.pi*d))

sympy.pprint(gamma_V)

cL = 0.6818817+0.6129331*d-0.1838439*d*d-0.1156844*d*d*d
cG = 0.3246017-0.6182531*d+0.1768139*d*d+0.1210944*d*d*d

Profile = cL*1./sympy.pi*gamma_V/((nu-nu_0)*(nu-nu_0)+gamma_V*gamma_V)+ \
          cG*sympy.sqrt(Log2)/(sympy.sqrt(sympy.pi)*gamma_L)*sympy.exp(-Log2*(nu-nu_0)**2/(gamma_V*gamma_V))

d_Profile__d_gamma_D = sympy.diff(Profile,gamma_D)
d_Profile__d_gamma_L = sympy.diff(Profile,gamma_L)

sympy.pprint(d_Profile__d_gamma_D)
sympy.pprint(d_Profile__d_gamma_L)
print("Starting the symbolic simplification::")

Profile =sympy.simplify(Profile)

input("Waiting here for the symbolic simplification::")
sympy.pprint(Profile)
