import math
from random import random
from numpy.linalg import norm as get_norm
import numpy as np
from main import elastic_collision

HYDROGEN_SCATTERING_PERCENTAGE = 7.8 / 10.5
ABSORBANCE_PERCENTAGE = 0.010063 / 0.361663
# Unitless, initial velocity 1MeV, thermalized velocity threshold 1eV
THERMALIZED_VELOCITY_THRESHOLD = 1 / math.sqrt(10 ** 6)
# In cm-1
MACROSCOPIC_CROSS_SECTION = 0.361663


def get_distance_to_next_interaction(macroscopic_cross_section):
    return - math.log(random()) / macroscopic_cross_section


def get_atomic_mass_target():
    if random() < (HYDROGEN_SCATTERING_PERCENTAGE):
        return 1
    else:
        return 16


def is_outside_right(position):
    if position[0] > 30:
        return True
    return False


def is_outside_left(position):
    if position[0] < 0:
        return True
    return False


def is_absorbed():
    return random() < ABSORBANCE_PERCENTAGE


def is_thermalized(v_f):
    return get_norm(v_f) < THERMALIZED_VELOCITY_THRESHOLD

def calculate_percentage(lst):
    unique_elements = set(lst)  # Get unique elements in the list
    total = len(lst)  # Total number of elements
    percentages = {element: (lst.count(element) / total) * 100 for element in unique_elements}  # Calculate percentages
    sorted_percentages = {k: percentages[k] for k in sorted(percentages)}  # Sort by key
    return sorted_percentages


# Possible simulation results: ESCAPED_RIGHT, ESCAPED_LEFT, ABSORBED, THERMALIZED

def simulate(x, v):
    while True:
        distance = get_distance_to_next_interaction(MACROSCOPIC_CROSS_SECTION)
        x = x + ((distance / get_norm(v_0)) * v_0)
        if is_outside_right(x):
            return "ESCAPED_RIGHT"

        if is_outside_left(x):
            return "ESCAPED_LEFT"

        if is_absorbed():
            return "ABSORBED"

        v = elastic_collision(v, get_atomic_mass_target())

        if is_thermalized(v):
            return "THERMALIZED"


v_0 = np.array([1, 0, 0])
# Position in cm
x_0 = np.array([0, 0, 0])

results = [simulate(x_0, v_0) for i in range(100000)]

print(calculate_percentage(results))
