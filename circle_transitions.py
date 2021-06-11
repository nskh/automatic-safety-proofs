from utils import *

r = 4
circle_traj = x ** 2 + y ** 2 - r ** 2  # radius 4
plot_implicit(circle_traj)

# TODO(nishant): plot shape on trajectory at given points

print(slope_sym(circle_traj))

transitions = find_transitions(circle_traj, angles)
transitions
