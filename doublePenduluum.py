import numpy as np
import sympy as sym
from sympy.abc import t
# %matplotlib inline
import matplotlib.pyplot as plt
from graphics import *
import time

#######################

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


def animate_double_pend(theta_array, L1 = 1, L2 = 1 ,T=10):
    
    windowWidth = 700.0
    windowHeight = 700.0

    # Create the objects that we will manipulate
    window = GraphWin("Double_Penduluum", windowWidth, windowHeight)
    window.setBackground("white")

    # Record the (x, y) start point in the window
    start_point = Point( (0.50) * (windowWidth), (0.70) * (windowHeight) )
    
    radius = 25.0

    Mass_1_Point = Point( 0, 0 )
    Mass_2_Point = Point( 0, 0 )
    Mass_1_Circle = Circle(Mass_1_Point, radius)
    Mass_2_Circle = Circle(Mass_2_Point, radius)

    link_1 = Line(start_point, Point(150, 100) )
    link_2 = Line(start_point, Point(50, 400) )
    
    # Draw Link-1
    link_1.draw(window)

    # Draw Mass-1
    Mass_1_Circle.setFill("red")
    Mass_1_Circle.draw(window)

    # Draw Link-2
    link_2.draw(window)
        
    # Draw Mass-2
    Mass_2_Circle.setFill("blue")
    Mass_2_Circle.draw(window)
    
    for i in range(len(theta_array) ):
        # Compute the next (x,y) for each component
        
        # line.move(10, 40)
        # undraw and redraw
        

        # Will need to tune this
        time.sleep(1.0)


 




# Create symbols 
g, m1, m2, theta1, theta2, R1 = sym.symbols('g m1 m2 theta1 theta2 R1')
R2, lambd = sym.symbols('R2 \lambda')

# Create the functions of time 
theta1 = sym.Function('theta1')(t)
theta2 = sym.Function('theta2')(t)
phi = sym.Function('phi')(theta1, theta2)


theta1_dot = theta1.diff(t)
theta1_dot_dot = theta1.diff(t, t)

theta2_dot = theta2.diff(t)
theta2_dot_dot = theta2.diff(t,t)

# First, write out the position of the two masses 
# We will differentiate these to get the mass's velocities
x1 = R1 * sym.sin(theta1)

y1 = -1 * (R1 * sym.cos(theta1) )


x2 = x1 + (R2 * sym.sin(theta1 + theta2) )

y2 =  (y1 - R2 * sym.cos(theta1 + theta2) )

x1_dt = x1.diff(t)
y1_dt = y1.diff(t)

x2_dt = x2.diff(t)
y2_dt = y2.diff(t)

phi = ( sym.sqrt( ( (x2**2) + (y2**2)) ) - np.sqrt(2) )


# Kinetic Energy of Mass 1
KE1 = (0.5) * (m1) * ( (x1_dt**2) + (y1_dt**2) )


# Potential Energy of Mass 1
# Note how we defined the height
V1 = m1 * g * y1


# Kinetic Energy of Mass2 
KE2 = (0.5) * (m2) * ( (x2_dt**2) + (y2_dt**2) )

# Potential Energy of Mass 2
V2 = m2 * g * y2

Lagrangian = KE1 + KE2 - V1 - V2
# Compute the Euler-Lagrange Equations
EL1 = (KE1 + KE2 - V1 - V2).diff(theta1_dot, t) - ( (KE1 + KE2 - V1 - V2).diff(theta1) )

EL2 = (KE1 + KE2 - V1 - V2).diff(theta2_dot, t) - ( (KE1 + KE2 - V1 - V2).diff(theta2) )

# Compute the gradient of phi
phi_gradient_1 = phi.diff(theta1)
phi_gradient_2 = phi.diff(theta2)

EL1 = sym.Eq( EL1, (lambd * phi_gradient_1) )
EL2 = sym.Eq( EL2, (lambd * phi_gradient_2) )


# Further differentiate phi
phi_dt = phi.diff(t)
phi_dt_dt = phi.diff(t, t)

# Remember to make phi_dt_dt an equation!
phi_dt_dt = sym.Eq(0, phi_dt_dt)

#####TESTING#########
#phi = sym.Eq(phi, 0)
#display(phi)
#display(phi_dt )
######TESTING########

# Substitute dummy values
a, b, c, d = sym.symbols('a b c d') 

EL1 = EL1.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d } ) # substitution

EL2 = EL2.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d } ) # substitution

phi_dt_dt = phi_dt_dt.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d } )


matrix_eq = sym.Matrix( [EL1, EL2, phi_dt_dt] )


# Solve the matrix for theta1_dot_dot and theta2_dot_dot
q = sym.Matrix( [b, d, lambd ] )


matrix_sol = sym.solve( matrix_eq, q )


subs1 = matrix_sol[b].subs(  {g: 9.81, R1: 1, R2: 1, m1: 1, m2:1  } )
subs2 = matrix_sol[d].subs(  {g: 9.81, R1: 1, R2: 1, m1: 1, m2:1  } )  
#subs3 = matrix_sol[lambd].subs(  {g: 9.81, R1: 1, R2: 1, m1: 1, m2:1  } ) 

computeTheta1_dt_dt = sym.lambdify( [theta1, a, theta2, c], subs1) 

computeTheta2_dt_dt = sym.lambdify( [theta1, a, theta2, c], subs2) 

#computeLambda = sym.lambdify( [theta1, a, theta2, c], subs3)



# This is my dynamics equation
def dynamics(q):
   
  theta1_dt_dt = computeTheta1_dt_dt( q[0], q[1], q[2], q[3]  )    
  
  theta2_dt_dt = computeTheta2_dt_dt( q[0], q[1], q[2], q[3]  )  
  
  return np.array( [q[1], theta1_dt_dt, q[3], theta2_dt_dt]   )



tspan = [0, 10]
dt = 0.01

# Define the initial conditions
theta1 = -np.pi/2.0
theta1_dot = 0.0

theta2 = -np.pi/2.0
theta2_dot = 0.0

initial_conditions = np.array( [theta1, theta1_dot, theta2, theta2_dot] )

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
plt.xlim(tspan)

plt.title('Double Penduluum Subject to Constraint')
plt.xlabel('Time (s)')
plt.ylabel('Angles Value')
plt.legend([r'$theta_1(t)$',r'$theta_2(t)$'])
plt.grid(True)
#plt.show()

print("Starting graphics")
theta_array = np.array( [ xvec[0], xvec[2] ]  )
animate_double_pend(theta_array, 1, 1, 10)

while (True):
    pass


