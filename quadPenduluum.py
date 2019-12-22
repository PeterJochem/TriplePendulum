import numpy as np
import sympy as sym
from sympy.abc import t
import matplotlib.pyplot as plt
from graphics import *

"""
    This function takes in an initial condition x0 and a timestep dt,
    as well as a dynamical system f(x) that outputs a vector of the
    same dimension as x0. It outputs a vector x at the future time step.
    Implements the Runge-Kutta Integration scheme
    """
def integrate(f,x0,dt):
    k1=dt*f(x0)
    k2=dt*f(x0+k1/2.)
    k3=dt*f(x0+k2/2.)
    k4=dt*f(x0+k3)
    xnew=x0+(1/6.)*(k1+2.*k2+2.*k3+k4)
    return xnew

  """
    This function takes in an initial condition x0, a timestep dt,
    a time span tspan consisting of a list [min_time, max_time],
    as well as a dynamical system f(x) that outputs a vector of the
    same dimension as x0. Additionally, this includes a flag (default false)
    that allows one to supply an Euler intergation scheme instead of
    the given scheme. It outputs a full trajectory simulated
    over the time span of dimensions (xvec_size, time_vec_size).
    """
def simulate(f,x0,tspan,dt,euler_int=False):
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

# This method animates the quadruple penduluum
# theta_array is a 4 X N array of the angle values over time 
def animate_quad_pend(theta_array):
 
    windowWidth = 700.0
    windowHeight = 700.0
    
    # Lists for changing the display
    # Could try showing only the last to or so points
    # history_1 = []
    # history_2 = []

    # For plotting
    history_radius = 1.0

    # Create the objects that we will manipulate
    window = GraphWin("Quadruple_Penduluum", windowWidth, windowHeight)
    window.setBackground("black")

    # Record the (x, y) start point in the window
    start_point = Point( (0.50) * (windowWidth), (0.30) * (windowHeight) )

    radius = 10.0

    Mass_1_Point = Point( 0, 0 )
    Mass_2_Point = Point( 0, 0 )
    Mass_3_Point = Point( 0, 0 )
    Mass_4_Point = Point( 0, 0 )

    Mass_1_Circle = Circle(Mass_1_Point, radius)
    Mass_2_Circle = Circle(Mass_2_Point, radius)
    Mass_3_Circle = Circle(Mass_3_Point, radius)
    Mass_4_Circle = Circle(Mass_4_Point, radius)

    link_1 = Line(start_point, Mass_1_Circle.getCenter() )
    link_2 = Line(Mass_1_Circle.getCenter(), Mass_2_Circle.getCenter() )
    link_3 = Line(Mass_2_Circle.getCenter(), Mass_3_Circle.getCenter() )
    link_4 = Line(Mass_3_Circle.getCenter(), Mass_4_Circle.getCenter() )

    # Draw Link-1
    link_1.draw(window) 
    # Draw Mass-1
    Mass_1_Circle.setFill("white")
    Mass_1_Circle.draw(window)
    Mass_1_Circle.setOutline("gray")

    # Draw Link-2
    link_2.draw(window)
    # Draw Mass-2
    Mass_2_Circle.setFill("white")
    Mass_2_Circle.draw(window)
    Mass_2_Circle.setOutline("gray")

    # Draw Link-3
    link_3.draw(window)
    # Draw Mass-3
    Mass_3_Circle.setFill("white")
    Mass_3_Circle.draw(window)
    Mass_3_Circle.setOutline("gray")
    
    # Draw Link-4
    link_4.draw(window)
    # Draw Mass-3
    Mass_4_Circle.setFill("white")
    Mass_4_Circle.draw(window)
    Mass_4_Circle.setOutline("gray")

    linkLength = 70.0
    
    for i in range(len(theta_array[0] ) ):

        # Compute the next (x,y) for each component
        x1_new = int( round(np.sin( theta_array[0][i] ) * linkLength) ) + start_point.getX()
        y1_new = int( round(np.cos( theta_array[0][i] ) * linkLength) ) + start_point.getY()

        delta_x_1 = x1_new - Mass_1_Circle.getCenter().getX()
        delta_y_1 = y1_new - Mass_1_Circle.getCenter().getY()

        x2_new = (np.sin( theta_array[0][i] + theta_array[1][i]  ) * linkLength) + Mass_1_Circle.getCenter().getX()
        y2_new = (np.cos( theta_array[0][i] + theta_array[1][i] ) * linkLength) + Mass_1_Circle.getCenter().getY()

        delta_x_2 = x2_new - Mass_2_Circle.getCenter().getX()
        delta_y_2 = y2_new - Mass_2_Circle.getCenter().getY()
        
        x3_new = (np.sin( theta_array[0][i] + theta_array[1][i] + theta_array[2][i]  ) * linkLength) + Mass_2_Circle.getCenter().getX()
        y3_new = (np.cos( theta_array[0][i] + theta_array[1][i] + theta_array[2][i]  ) * linkLength) + Mass_2_Circle.getCenter().getY()
        
        delta_x_3 = x3_new - Mass_3_Circle.getCenter().getX()
        delta_y_3 = y3_new - Mass_3_Circle.getCenter().getY()
        
        x4_new = (np.sin( theta_array[0][i] + theta_array[1][i] + theta_array[2][i] + theta_array[3][i]  ) * linkLength) + Mass_3_Circle.getCenter().getX()
        y4_new = (np.cos( theta_array[0][i] + theta_array[1][i] + theta_array[2][i] + theta_array[3][i]  ) * linkLength) + Mass_3_Circle.getCenter().getY()

        delta_x_4 = x4_new - Mass_4_Circle.getCenter().getX()
        delta_y_4 = y4_new - Mass_4_Circle.getCenter().getY()
        
        # Undraw the first link 
        link_1.undraw()
        # Draw the first link
        link_1 =  Line(start_point, Mass_1_Circle.getCenter() )
        link_1.draw(window)
        link_1.setFill("gray")

        # Undraw the first link
        link_2.undraw()
        # Draw the second link  
        link_2 = Line(Mass_1_Circle.getCenter(), Mass_2_Circle.getCenter() )
        link_2.draw(window)
        link_2.setFill("gray")
        
        # Undraw the first link
        link_3.undraw()
        # Draw the second link
        link_3 = Line(Mass_2_Circle.getCenter(), Mass_3_Circle.getCenter() )
        link_3.draw(window)
        link_3.setFill("gray")

        link_4.undraw()
        # Draw the second link
        link_4 = Line(Mass_3_Circle.getCenter(), Mass_4_Circle.getCenter() )
        link_4.draw(window)
        link_4.setFill("gray")

        # Draw the next frame
        Mass_1_Circle.move(delta_x_1, delta_y_1)
        Mass_2_Circle.move(delta_x_2, delta_y_2)
        Mass_3_Circle.move(delta_x_3, delta_y_3)
        Mass_4_Circle.move(delta_x_4, delta_y_4)

        # Add the plot to the history
        history_1 = Circle( Point( Mass_1_Circle.getCenter().getX(), Mass_1_Circle.getCenter().getY() ) , history_radius)
        history_2 = Circle( Point( Mass_2_Circle.getCenter().getX(), Mass_2_Circle.getCenter().getY() ) , history_radius)
        history_3 = Circle( Point( Mass_3_Circle.getCenter().getX(), Mass_3_Circle.getCenter().getY() ) , history_radius)
        history_4 = Circle( Point( Mass_4_Circle.getCenter().getX(), Mass_4_Circle.getCenter().getY() ) , history_radius)
        # Draw the history points
        # Delay this? 
        history_1.draw(window)
        history_2.draw(window)
        history_3.draw(window)
        history_1.setFill("blue")
        history_2.setFill("red")
        history_3.setFill("orange")
        history_4.setFill("green")
        
        # Re-draw the masses 
        # Draw Mass-1
        Mass_1_Circle.undraw()
        Mass_1_Circle.draw(window)

        # Draw Mass-2
        Mass_2_Circle.undraw()
        Mass_2_Circle.draw(window)

        # Draw Mass-3
        Mass_3_Circle.undraw()
        Mass_3_Circle.draw(window)
        
        # Draw Mass-4
        Mass_4_Circle.undraw()
        Mass_4_Circle.draw(window)

        # Will need to tune this
        #           0.005
        time.sleep(0.00005)                                                                       


# Create symbols 
g, m, R, theta1, theta2, theta3, theta4 = sym.symbols('g m R theta1 theta2 theta3 theta4')

# Create the functions of time 
theta1 = sym.Function('theta1')(t)
theta2 = sym.Function('theta2')(t)
theta3 = sym.Function('theta3')(t)
theta4 = sym.Function('theta4')(t)

theta1_dot = theta1.diff(t)
theta1_dot_dot = theta1.diff(t, t)

theta2_dot = theta2.diff(t)
theta2_dot_dot = theta2.diff(t, t)

theta3_dot = theta3.diff(t)
theta3_dot_dot = theta3.diff(t, t)

theta4_dot = theta4.diff(t)
theta4_dot_dot = theta4.diff(t, t)

# First, write out the position of the two masses 
# We will differentiate these to get the mass's velocities
x1 = R * sym.sin(theta1)
y1 = -1 * (R * sym.cos(theta1) )

x2 = x1 + (R * sym.sin(theta1 + theta2) )
y2 = y1 - (R * sym.cos(theta1 + theta2) )

x3 = x2 + (R * sym.sin(theta1 + theta2 + theta3) )
y3 = y2 - (R * sym.cos(theta1 + theta2 + theta3) )

x4 = x3 + (R * sym.sin(theta1 + theta2 + theta3 + theta4) )
y4 = y3 - (R * sym.cos(theta1 + theta2 + theta3 + theta4) )

x1_dt = x1.diff(t)
y1_dt = y1.diff(t)

x2_dt = x2.diff(t)
y2_dt = y2.diff(t)

x3_dt = x3.diff(t)
y3_dt = y3.diff(t)

x4_dt = x4.diff(t)
y4_dt = y4.diff(t)

# Kinetic Energy of Mass 1
KE1 = (0.5) * (m) * ( (x1_dt**2) + (y1_dt**2) )
# Potential Energy of Mass 1
# Note how we defined the height
V1 = m * g * y1

# Kinetic Energy of Mass 2  
KE2 = (0.5) * (m) * ( (x2_dt**2) + (y2_dt**2) )
# Potential Energy of Mass 2
V2 = m * g * y2

# Kinetic Energy of Mass 3 
KE3 = (0.5) * (m) * ( (x3_dt**2) + (y3_dt**2) )
# Potential Energy of Mass 3
V3 = m * g * y3

# Kinetic Energy of Mass 4
KE4 = (0.5) * (m) * ( (x4_dt**2) + (y4_dt**2) )
# Potential Energy of Mass 4
V4 = m * g * y4

Lagrangian = (KE1 + KE2 + KE3 + KE4) - (V1 + V2 + V3 + V4)
# Compute the Euler-Lagrange Equations
EL1 = Lagrangian.diff(theta1_dot, t) - ( Lagrangian.diff(theta1) )

EL2 = Lagrangian.diff(theta2_dot, t) - ( Lagrangian.diff(theta2) )

EL3 = Lagrangian.diff(theta3_dot, t) - ( Lagrangian.diff(theta3) )

EL4 = Lagrangian.diff(theta4_dot, t) - ( Lagrangian.diff(theta4) )

# Call sym.simplify
EL1 = sym.simplify(EL1)
EL2 = sym.simplify(EL2)
EL3 = sym.simplify(EL3)
EL4 = sym.simplify(EL4)

EL1 = sym.Eq( EL1, 0.0 )
EL2 = sym.Eq( EL2, 0.0 )
EL3 = sym.Eq( EL3, 0.0 )
EL4 = sym.Eq( EL4, 0.0 )

# Substitute dummy values
a, b, c, d, e, f, g, h = sym.symbols('a b c d e f g h') 

EL1 = EL1.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f, theta4.diff(t): g, theta4.diff(t, t): h} ) 

EL2 = EL2.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f, theta4.diff(t): g, theta4.diff(t, t): h} )

EL3 = EL3.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f, theta4.diff(t): g, theta4.diff(t, t): h} ) 

EL4 = EL4.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f, theta4.diff(t): g, theta4.diff(t, t): h} )

matrix_eq = sym.Matrix( [EL1, EL2, EL3, EL4] )

# Solve the matrix for theta1_dot_dot and theta2_dot_dot
desiredVars = sym.Matrix( [b, d, f, h] )

matrix_sol = sym.solve( matrix_eq, desiredVars )

print("The matrix solution is")
print(matrix_sol)

subs1 =  matrix_sol[b].subs(  {g: 9.81, R: 1, m: 1} )
subs2 =  matrix_sol[d].subs(  {g: 9.81, R: 1, m: 1} )  
subs3 =  matrix_sol[f].subs(  {g: 9.81, R: 1, m: 1} )
subs4 =  matrix_sol[h].subs(  {g: 9.81, R: 1, m: 1} )

computeTheta1_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e, theta4, g], subs1) 

computeTheta2_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e, theta4, g], subs2) 

computeTheta3_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e, theta4, g], subs3) 

computeTheta4_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e, theta4, g], subs4)

# This is my dynamics equation
# It computes the new accelerations 
def dynamics(q):
   
  theta1_dt_dt = computeTheta1_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7] )    
  
  theta2_dt_dt = computeTheta2_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7] )  
  
  theta3_dt_dt = computeTheta3_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7] )
    
  theta4_dt_dt = computeTheta4_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7] )

  return np.array( [q[1], theta1_dt_dt, q[3], theta2_dt_dt, q[5], theta3_dt_dt, q[7], theta4_dt_dt] )



# Define the initial conditions
tspan = [0, 100]
dt = 0.01

theta1 = (np.pi / 2.0) * (4.0 / 3.0)
theta1_dot = 0.0

theta2 = 0.01
theta2_dot = 0.0

theta3 = 0.1
theta3_dot = 0.0

theta4 = 0.01
theta4_dot = 0.0

initial_conditions = np.array( [theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, theta4, theta4_dot] )

N = int( (max(tspan) - min(tspan) ) / dt )
tvec = np.linspace(min(tspan), max(tspan), N)

##################
# xvec is an array of each angle's value over time
xvec = simulate(dynamics, initial_conditions, tspan, dt) 
##############

plt.figure(dpi = 110, facecolor = 'w')
# Plot each of the angle positions over time
plt.plot( tvec, xvec[0] )
plt.plot( tvec, xvec[2] )
plt.plot( tvec, xvec[4] )
plt.plot( tvec, xvec[6] )
plt.xlim(tspan)

plt.title('Quadruple Penduluum')
plt.xlabel('Time (s)')
plt.ylabel('Angles Value')
plt.legend([r'$theta_1(t)$',r'$theta_2(t)$', r'$theta_3(t)$', r'$theta_4(t)$'])
plt.grid(True)
# plt.show()

g = input("Press Enter to animate the quadruple penduluum")

theta_array = np.array( [ xvec[0], xvec[2], xvec[4], xvec[6] ] )
animate_quad_pend(theta_array)



