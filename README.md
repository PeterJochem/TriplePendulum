How to run my code
To run the double penduluum simulation
```python3 doublePenduluum.py``` (Runtim prior to animation is about a few seconds)
To run the triple penduluum simulation (Runtime prior to animation is about 1 minute)
```python3 TriplePenduluum.py```
To run the quadruple penduluum simulation (Runtime prior to animation is many hours)
```python3 quadPenduluum```
To see multiple triple penduluums with slightly diffrent initial conditions run ```python3 chaos.py``` (Runtime prior to animation is about 1 minute)


Results\n
The triple penduluum is extremly sensitive to perturbations!! Don't beleive me? Run my code and see how quickly slight diffrences in the initial conditions leads to very diffrent trajectories.

Add T_0_Image

T_1_Second_Image

T_2_Second_Image

Insert Images  

A really awesome surprise was how cool the trajectory plotted over time looked! Heres a pic. This could realy make for a cool screensaver

![Triple Penduluum](https://github.com/PeterJochem/TriplePendulum/blob/master/Trajectory.png)

FIX Me: Add video to the repo

Chaos



Improvements\n
I used sympy to solve the equations of motion symbolically. This is fine for the single, double, and triple penduluum but becomes an unreasonable approach when you go beyond three penduluums. I should try solving numerically and see how the performace compares. It should be dramatically faster.
