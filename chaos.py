import numpy as np
import sympy as sym
from sympy.abc import t
import matplotlib.pyplot as plt
from graphics import *


# This class runs multiple penduluums to see how chaos ensues!
# A slight perturbation in the inital conditons leads to dramaticly diffrent 
# paths through space
class chaos:
    
    # Describe the input parameters
    def __init__(self, numPends, error):
        
        # This is the graphics window we will wrtie to
        # self.window =
        
        self.computeTheta1_dt_dt = None
        self.computeTheta2_dt_dt = None
        self.computeTheta3_dt_dt = None
        
        # Derive the equations of motion
        self.create_Equations_of_Motion()
        
        self.numPends = numPends
        self.error = error
        # This is each penduluums array of positions over time
        self.all_penduluum_positions = [ ]         
        
        for i in range(self.numPends):
                        
            # Create the next initial conditions
            theta1 = (np.pi / 2.0) * (4.0 / 3.0) + (i * error)
            theta1_dot = 0.0

            theta2 = -1 * (np.pi / 8.0)
            theta2_dot = 0.0

            theta3 = np.pi / 2.0
            theta3_dot = 0

            initial_conditions = np.array( [theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot] )

            pend_next = self.generate_penduluum_positions( initial_conditions, 0.01 )
            self.all_penduluum_positions.append(pend_next) 
            
        # Reorganizes the all_penduluum_positions array to make graphics code simpler
        self.reshapeArray()
        # theta_array =  np.array( [ self.all_penduluum_positions[0][0], self.all_penduluum_positions[0][2], self.all_penduluum_positions[0][4] ] )
        self.animate_multiple_pend( self.all_penduluum_positions )
         

    """
    This function takes in an initial condition x0 and a timestep dt,
    as well as a dynamical system f(x) that outputs a vector of the
    same dimension as x0. It outputs a vector x at the future time step.
    """
    def integrate(self, f, x0, dt):
        k1 = dt * f(x0)
        k2 = dt * f(x0 + k1/2.)
        k3 = dt * f(x0 + k2/2.)
        k4 = dt * f(x0 + k3)
        xnew = x0 + (1/6.) * (k1 + 2. * k2 + 2.* k3 + k4)
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
    def simulate(self, f, x0, tspan, dt, euler_int = False):
        N = int((max(tspan)-min(tspan))/dt)
        x = np.copy(x0) 
        tvec = np.linspace(min(tspan),max(tspan),N)
        xtraj = np.zeros((len(x0),N))
        for i in range(N):
            if euler_int:
                xtraj[:,i] = self.euler(f,x,dt)
            else:
                xtraj[:,i] = self.integrate(f,x,dt)
            x = np.copy(xtraj[:,i])
        return xtraj   
    
    # This method reorganizes the all_penduluum_positions array
    # to make writing the graphics code easier
    def reshapeArray(self):
        
        old_array = self.all_penduluum_positions.copy()
        
        new_array = [ ]

        for i in range( self.numPends ):
            # theta_array =  np.array( [ self.all_penduluum_positions[0][0], self.all_penduluum_positions[0][2], self.all_penduluum_positions[0][4] ] )
            next_entry = np.array( [ self.all_penduluum_positions[i][0], self.all_penduluum_positions[i][2], self.all_penduluum_positions[i][4] ] ) 
            new_array.append( next_entry )
            
        self.all_penduluum_positions = new_array 


    # Describe here
    def animate_multiple_pend(self, theta_array):
        
        windowWidth = 700.0
        windowHeight = 700.0

        # For plotting
        history_radius = 1.4

        # Create the objects that we will manipulate
        window = GraphWin("Chaos", windowWidth, windowHeight)
        window.setBackground("black")

        # Record the (x, y) start point in the window
        start_point = Point( (0.50) * (windowWidth), (0.30) * (windowHeight) )
        radius = 10.0
        
        Mass_1_Point = []
        Mass_2_Point = []
        Mass_3_Point = []
        
        Mass_1_Circle = []
        Mass_2_Circle = []
        Mass_3_Circle = []

        link_1 = []
        link_2 = []
        link_3 = []

        mass_colors = [ "white", "red", "orange", "blue" ]

        for i in range(self.numPends):
            Mass_1_Point.append( Point( 0, 0 ) )
            Mass_2_Point.append( Point( 0, 0 ) )
            Mass_3_Point.append( Point( 0, 0 ) )
    
            Mass_1_Circle.append( Circle(Mass_1_Point[i], radius) )
            Mass_2_Circle.append( Circle(Mass_2_Point[i], radius) )
            Mass_3_Circle.append( Circle(Mass_3_Point[i], radius) )

            link_1.append( Line(start_point, Mass_1_Circle[i].getCenter() ) )
            link_2.append( Line(Mass_1_Circle[i].getCenter(), Mass_2_Circle[i].getCenter() ) )
            link_3.append( Line(Mass_2_Circle[i].getCenter(), Mass_3_Circle[i].getCenter() ) )

            # Draw Link-1
            link_1[i].draw(window)

            # Draw Mass-1
            Mass_1_Circle[i].setFill(mass_colors[i])
            Mass_1_Circle[i].draw(window)
            Mass_1_Circle[i].setOutline("gray")

            # Draw Link-2
            link_2[i].draw(window)

            # Draw Mass-2
            Mass_2_Circle[i].setFill(mass_colors[i])
            Mass_2_Circle[i].draw(window)
            Mass_2_Circle[i].setOutline("gray")

            # Draw Link-3
            link_3[i].draw(window)

            # Draw Mass-3
            Mass_3_Circle[i].setFill(mass_colors[i])
            Mass_3_Circle[i].draw(window)
            Mass_3_Circle[i].setOutline("gray")

        linkLength = 100.0
        
        # theta_array is a list of 2-D numpy arrays 
        for i in range(len(theta_array[0][0] ) ):                
            for j in range(self.numPends):
                # Compute the next (x,y) for each component
                x1_new = int( round(np.sin( theta_array[j][0][i] ) * linkLength) ) + start_point.getX()
                y1_new = int( round(np.cos( theta_array[j][0][i] ) * linkLength) ) + start_point.getY()

                delta_x_1 = x1_new - Mass_1_Circle[j].getCenter().getX()
                delta_y_1 = y1_new - Mass_1_Circle[j].getCenter().getY()

                x2_new = (np.sin( theta_array[j][0][i] + theta_array[j][1][i]  ) * linkLength) + Mass_1_Circle[j].getCenter().getX()
                y2_new = (np.cos( theta_array[j][0][i] + theta_array[j][1][i] ) * linkLength) + Mass_1_Circle[j].getCenter().getY()

                delta_x_2 = x2_new - Mass_2_Circle[j].getCenter().getX()
                delta_y_2 = y2_new - Mass_2_Circle[j].getCenter().getY()
        
                x3_new = (np.sin( theta_array[j][0][i] + theta_array[j][1][i] + theta_array[j][2][i]  ) * linkLength) + Mass_2_Circle[j].getCenter().getX()
                y3_new = (np.cos( theta_array[j][0][i] + theta_array[j][1][i] + theta_array[j][2][i]  ) * linkLength) + Mass_2_Circle[j].getCenter().getY()
        
                delta_x_3 = x3_new - Mass_3_Circle[j].getCenter().getX()
                delta_y_3 = y3_new - Mass_3_Circle[j].getCenter().getY()

                # Undraw the first link 
                link_1[j].undraw()
                # Draw the first link
                link_1[j] = Line(start_point, Mass_1_Circle[j].getCenter() )
                link_1[j].draw(window)
                link_1[j].setFill("gray")

                # Undraw the first link
                link_2[j].undraw()
                # Draw the second link  
                link_2[j] = Line(Mass_1_Circle[j].getCenter(), Mass_2_Circle[j].getCenter() )
                link_2[j].draw(window)
                link_2[j].setFill("gray")
        
                # Undraw the first link
                link_3[j].undraw()
                # Draw the second link
                link_3[j] = Line(Mass_2_Circle[j].getCenter(), Mass_3_Circle[j].getCenter() )
                link_3[j].draw(window)
                link_3[j].setFill("gray")

                # Draw the next frame
                Mass_1_Circle[j].move(delta_x_1, delta_y_1)
                Mass_2_Circle[j].move(delta_x_2, delta_y_2)
                Mass_3_Circle[j].move(delta_x_3, delta_y_3)

                # Add the plot to the history
                history_1 = Circle( Point( Mass_1_Circle[j].getCenter().getX(), Mass_1_Circle[j].getCenter().getY() ) , history_radius)
                history_2 = Circle( Point( Mass_2_Circle[j].getCenter().getX(), Mass_2_Circle[j].getCenter().getY() ) , history_radius)
                history_3 = Circle( Point( Mass_3_Circle[j].getCenter().getX(), Mass_3_Circle[j].getCenter().getY() ) , history_radius)
                # Draw the history points
                # Delay this? 
                #history_1.draw(window)
                #history_2.draw(window)
                history_3.draw(window)
                history_3.setFill(mass_colors[j]) 
                #history_2.setFill("red")
                #history_3.setFill("orange")
        
                # Re-draw the masses 
                # Draw Mass-1
                #Mass_1_Circle[j].undraw()
                #Mass_1_Circle[j].draw(window)

                # Draw Mass-2
                #Mass_2_Circle[j].undraw()
                #Mass_2_Circle[j].draw(window)

                # Draw Mass-3
                #Mass_3_Circle[j].undraw()
                #Mass_3_Circle[j].draw(window)

            # Will need to tune this
            #           0.005
            time.sleep(0.05)                                                                       

    def create_Equations_of_Motion(self):
        # Create symbols 
        g, m, R, theta1, theta2, theta3 = sym.symbols('g m R theta1 theta2 theta3')

        # Create the functions of time 
        theta1 = sym.Function('theta1')(t)
        theta2 = sym.Function('theta2')(t)
        theta3 = sym.Function('theta3')(t)

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
        y2 = (y1 - R * sym.cos(theta1 + theta2) )

        x3 = x2 + (R * sym.sin(theta1 + theta2 + theta3) )
        y3 = y2 - (R * sym.cos(theta1 + theta2 + theta3) )

        x1_dt = x1.diff(t)
        y1_dt = y1.diff(t)

        x2_dt = x2.diff(t)
        y2_dt = y2.diff(t)

        x3_dt = x3.diff(t)
        y3_dt = y3.diff(t)

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

        EL1 = sym.Eq( EL1, 0.0 )
        EL2 = sym.Eq( EL2, 0.0 )
        EL3 = sym.Eq( EL3, 0.0 )

        # Substitute dummy values
        a, b, c, d, e, f = sym.symbols('a b c d e f') 

        EL1 = EL1.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f } )
        EL2 = EL2.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f } )
        EL3 = EL3.subs( { theta1.diff(t): a, theta1.diff(t,t): b, theta2.diff(t): c, theta2.diff(t,t): d, theta3.diff(t): e, theta3.diff(t,t): f } )
        matrix_eq = sym.Matrix( [EL1, EL2, EL3] )

        # Solve the matrix for theta1_dot_dot and theta2_dot_dot
        desiredVars = sym.Matrix( [b, d, f] )
        matrix_sol = sym.solve( matrix_eq, desiredVars )

        subs1 = matrix_sol[b].subs( {g: 9.81, R: 1, m: 1} )
        subs2 = matrix_sol[d].subs( {g: 9.81, R: 1, m: 1} )  
        subs3 = matrix_sol[f].subs( {g: 9.81, R: 1, m: 1} )

        # Keep the lambdified expressions in the class structure to use later 
        self.computeTheta1_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e], subs1) 
        self.computeTheta2_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e], subs2) 
        self.computeTheta3_dt_dt = sym.lambdify( [theta1, a, theta2, c, theta3, e], subs3) 

    # This is my dynamics equation
    def dynamics(self, q):
   
        theta1_dt_dt = self.computeTheta1_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5] )    
  
        theta2_dt_dt = self.computeTheta2_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5] )  
 
        theta3_dt_dt = self.computeTheta3_dt_dt( q[0], q[1], q[2], q[3], q[4], q[5] )
  
        return np.array( [q[1], theta1_dt_dt, q[3], theta2_dt_dt, q[5], theta3_dt_dt]   )

    # Take
    def generate_penduluum_positions( self, initial_conditions, dt ):

        # Define the initial conditions
        tspan = [0, 100]
        
        N = int( (max(tspan) - min(tspan) ) / dt )
        tvec = np.linspace(min(tspan), max(tspan), N)

        # Here we simulate
        xvec = self.simulate( self.dynamics, initial_conditions, tspan, dt) 
         
        return xvec


myChaos = chaos( 4, 0.01 )


# Create chaos object 
# theta_array = np.array( [ xvec[0], xvec[2], xvec[4] ]  )
# animate_triple_pend(theta_array)




