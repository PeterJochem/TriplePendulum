# WOOOO
import numpy as np
import sympy as sym
from sympy.abc import t
#%matplotlib inline
import matplotlib.pyplot as plt

#######################
# Custom latex printing
def custom_latex_printer(exp,**options):
    from google.colab.output._publish import javascript
    url = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.3/latest.js?config=TeX-AMS_HTML"
    javascript(url=url)
    return sym.printing.latex(exp,**options)
sym.init_printing(use_latex="mathjax",latex_printer=custom_latex_printer)

####################
# Simulation helpers
def integrate(f,x0,dt):
    """
    This function takes in an initial condition x0 and a timestep dt,
    as well as a dynamical system f(x) that outputs a vector of the
    same dimension as x0. It outputs a vector x at the future time step.
    """
    k1=dt*f(x0)
    k2=dt*f(x0+k1/2.)
    k3=dt*f(x0+k2/2.)
    k4=dt*f(x0+k3)
    xnew=x0+(1/6.)*(k1+2.*k2+2.*k3+k4)
    return xnew

def simulate(f,x0,tspan,dt,euler_int=False):
    """
    This function takes in an initial condition x0, a timestep dt,
    a time span tspan consisting of a list [min_time, max_time],
    as well as a dynamical system f(x) that outputs a vector of the
    same dimension as x0. Additionally, this includes a flag (default false)
    that allows one to supply an Euler intergation scheme instead of 
    the given scheme. It outputs a full trajectory simulated
    over the time span of dimensions (xvec_size, time_vec_size).
    """
    N = int((max(tspan)-min(tspan))/dt)
    x = np.copy(x0)
    tvec = np.linspace(min(tspan),max(tspan),N)
    xtraj = np.zeros((len(x0),N))
    for i in range(N):
        if euler_int:
            xtraj[:,i]=euler(f,x,dt)
        else:
            xtraj[:,i]=integrate(f,x,dt)
        x = np.copy(xtraj[:,i])
    return xtraj   


def animate_triple_pend(theta_array,L1=1,L2=1,L3=1,T=10):
    
    ################################
    # Imports required for animation.
    from plotly.offline import init_notebook_mode, iplot
    from IPython.display import display, HTML

    #######################
    # Browser configuration.
    def configure_plotly_browser_state():
      import IPython
      display(IPython.core.display.HTML('''
            <script src="/static/components/requirejs/require.js"></script>
            <script>
              requirejs.config({
                paths: {
                  base: '/static/base',
                  plotly: 'https://cdn.plot.ly/plotly-1.5.1.min.js?noext',
                },
              });
            </script>
            '''))
    configure_plotly_browser_state()
    init_notebook_mode(connected=False)

    ###############################################
    # Getting data from pendulum angle trajectories.
    xx1 = L1*np.sin(theta_array[0])
    yy1 = -L1*np.cos(theta_array[0])

    xx2 = xx1 + L2 * np.sin( theta_array[0] + theta_array[1] )
    yy2 = yy1 - L2 * np.cos( theta_array[0] + theta_array[1] )
    
    xx3 = xx2 +  L3 * np.sin( theta_array[0] + theta_array[1] + theta_array[2] )
    yy3 = yy2 - L3 * np.cos( theta_array[0] + theta_array[1] + theta_array[2]  )


    N = len(theta_array[0]) # Need this for specifying length of simulation

    ####################################
    # Using these to specify axis limits.
    xm=np.min(xx1)-0.5
    xM=np.max(xx1)+0.5
    ym=np.min(yy1)-2.5
    yM=np.max(yy1)+1.5


# Create symbols 
g, m, R, theta1, theta2, theta3 = sym.symbols('g m R theta1 theta2 theta3')
lambd = sym.symbols('\lambda')

# Trying subbing in/ just using 1's for all R and m

# Create the functions of time 
theta1 = sym.Function('theta1')(t)
theta2 = sym.Function('theta2')(t)
theta3 = sym.Function('theta3')(t)

phi = sym.Function('phi')(theta1, theta2, theta3)


theta1_dot = theta1.diff(t)
theta1_dot_dot = theta1.diff(t, t)

theta2_dot = theta2.diff(t)
theta2_dot_dot = theta2.diff(t,t)

theta3_dot = theta3.diff(t)
theta3_dot_dot = theta3.diff(t,t)

# First, write out the position of the two masses 
# We will differentiate these to get the mass's velocities
x1 = R * sym.sin(theta1)

y1 = -1 * (R * sym.cos(theta1) )


x2 = x1 + (R * sym.sin(theta1 + theta2) )
y2 =  (y1 - R * sym.cos(theta1 + theta2) )

x3 = x2 + (R * sym.sin(theta1 + theta2 + theta3) )

y3 = y2 - (R * sym.cos(theta1 + theta2 + theta3) )


x1_dt = x1.diff(t)
y1_dt = y1.diff(t)

x2_dt = x2.diff(t)
y2_dt = y2.diff(t)

x3_dt = x3.diff(t)
y3_dt = y3.diff(t)


phi =  (x3**2) + (y3**2) - 5


# Kinetic Energy of Mass 1
KE1 = (0.5) * (m) * ( (x1_dt**2) + (y1_dt**2) )


# Potential Energy of Mass 1
# Note how we defined the height
V1 = m * g * y1

# Kinetic Energy of Mass2 
KE2 = (0.5) * (m) * ( (x2_dt**2) + (y2_dt**2) )

# Potential Energy of Mass 2
V2 = m * g * y2

# Kinetic Energy of Mass3 
KE3 = (0.5) * (m) * (  (x3_dt**2) + (y3_dt**2) )

# Potential Energy of Mass 3
V3 = m * g * y3


Lagrangian = KE1 + KE2 + KE3 - V1 - V2 - V3
# Compute the Euler-Lagrange Equations
EL1 = (KE1 + KE2 + KE3 - V1 - V2 - V3).diff(theta1_dot, t) - ( (KE1 + KE2 + KE3 - V1 - V2 - V3).diff(theta1) )

EL2 = (KE1 + KE2 + KE3 - V1 - V2 - V3).diff(theta2_dot, t) - ( (KE1 + KE2 + KE3 - V1 - V2 - V3).diff(theta2) )

EL3 = (KE1 + KE2 + KE3 - V1 - V2 - V3).diff(theta3_dot, t) - ( (KE1 + KE2 + KE3 - V1 - V2 - V3).diff(theta3) )

# Call sym.simplify
EL1 = sym.simplify(EL1)
EL2 = sym.simplify(EL2)
EL3 = sym.simplify(EL3)


# Compute the gradient of phi
phi_gradient_1 = phi.diff(theta1)
phi_gradient_2 = phi.diff(theta2)
phi_gradient_3 = phi.diff(theta3)

# More simplifying
phi_gradient_1 = sym.simplify(phi_gradient_1)
phi_gradient_2 = sym.simplify(phi_gradient_2)
phi_gradient_3 = sym.simplify(phi_gradient_3)


EL1 = sym.Eq( EL1, (lambd * phi_gradient_1) )
EL2 = sym.Eq( EL2, (lambd * phi_gradient_2) )
EL3 = sym.Eq( EL3, (lambd * phi_gradient_3) )

# Further differentiate phi
phi_dt = phi.diff(t)
phi_dt_dt = phi.diff(t, t)

phi_dt_dt = sym.simplify(phi_dt_dt)

# Remember to make phi_dt_dt an equation!
phi_dt_dt = sym.Eq(0, phi_dt_dt)

#####TESTING#########
#phi = sym.Eq(phi, 0)
#display(phi)
#display(phi_dt )
######TESTING########

# Substitute dummy values
a, b, c, d, e, f = sym.symbols('a b c d e f') 

EL1 = EL1.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f } ) # substitution

EL2 = EL2.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f } ) # substitution

EL3 = EL3.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f } ) # substitution


phi_dt_dt = phi_dt_dt.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f } )


matrix_eq = sym.Matrix( [EL1, EL2, EL3, phi_dt_dt] )


# Solve the matrix for theta1_dot_dot and theta2_dot_dot
desiredVars = sym.Matrix( [b, d, f, lambd ] )

#display(matrix_eq)

matrix_sol = sym.solve( matrix_eq, desiredVars )

subs1 = matrix_sol[b].subs(  {g: 9.81, R: 1, m: 1} )
subs2 = matrix_sol[d].subs( {g: 9.81, R: 1, m: 1}   )  
subs3 =  matrix_sol[f].subs(  {g: 9.81, R: 1, m: 1}  )
#subs4 = matrix_sol[lambd].subs(  {g: 9.81, R1: 1, R2: 1, m1: 1, m2:1  } ) 

computeTheta1_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e], subs1) 

computeTheta2_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e], subs2) 

computeTheta3_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e], subs3) 

#computeLambda = sym.lambdify( [theta1, a, theta2, c], subs3)




# This is my dynamics equation
def dynamics(q):
   
  theta1_dt_dt = computeTheta1_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5] )    
  
  theta2_dt_dt = computeTheta2_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5] )  
  
  theta3_dt_dt = computeTheta3_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5] )
  
  return np.array( [q[1], theta1_dt_dt, q[3], theta2_dt_dt, q[5], theta3_dt_dt]   )



tspan = [0, 10]
dt = 0.01

# Define the initial conditions
theta1 = -np.pi/2.0
theta1_dot = 0.0

theta2 = -np.pi/2.0
theta2_dot = 0.0

theta3 = np.pi / 2.0
theta3_dot = 0

initial_conditions = np.array( [theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot] )

N = int( (max(tspan) - min(tspan) ) / dt )
tvec = np.linspace(min(tspan), max(tspan), N)

##################
# Here we simulate
xvec = simulate(dynamics, initial_conditions, tspan, dt) 
##############
# Here we plot
plt.figure(dpi=110,facecolor='w')

#plt.plot(tvec, xvec[0])
#plt.plot(tvec, xvec[1])
plt.plot(tvec, xvec[0] )
plt.plot(tvec, xvec[2] )
plt.plot( tvec, xvec[4] )
plt.xlim(tspan)

plt.title('Double Penduluum Subject to Constraint')
plt.xlabel('Time (s)')
plt.ylabel('Angles Value')
plt.legend([r'$theta_1(t)$',r'$theta_2(t)$', r'$theta_3(t)$'])
plt.grid(True)
plt.show()

#theta_array = np.array( [ xvec[0], xvec[2], xvec[4] ]  )
#animate_triple_pend(theta_array,L1=1,L2=1,L3=1,T=10)

theta_array = np.array( [ xvec[0], xvec[2], xvec[4] ]  )
animate_triple_pend(theta_array,L1=1,L2=1,L3=1,T=10)




