from sympy import *
from typing import Tuple, List, Dict, Set


class TransitionPoint:
    # TODO(nishant): add documentation here
    # May not always represent a transition point (sometimes is piecewise bound)
    def __init__(self, point: Point, angle, function):
        self.point = point
        self.angle = angle
        self.function = function
        self.is_bound = not angle  # used to see if this is a piecewise bound or not

    def __repr__(self):
        if self.is_bound:
            # return f"Boundary({self.point.x}, {self.point.y}) [func {self.function}]"
            return f"Boundary({self.point.x}, {self.point.y})"
        else:
            # return f"Transition({self.point.x}, {self.point.y}) [angle {self.angle}, func {self.function}]"
            return f"Transition({self.point.x}, {self.point.y})"


class TransitionPointSet:
    def __init__(self):
        self.objects: Dict = {}
        self.unique_points: set = set()

    def add(self, tp: TransitionPoint):
        if tp.point in self.unique_points:
            return
        else:
            self.objects[tp.point] = tp
            self.unique_points.add(tp.point)

    def update(self, tps: set):
        for tp in tps:
            self.add(tp)

    # TODO: make sortable
