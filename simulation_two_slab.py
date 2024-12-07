import math
from random import random

from numpy.linalg import norm as get_norm

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
    if random() < HYDROGEN_SCATTERING_PERCENTAGE:
        return 1
    else:
        return 16


def is_outside_right(position):
    return position[0] > 30


def is_outside_left(position):
    return position[0] < 0


def is_absorbed():
    return random() < ABSORBANCE_PERCENTAGE


def is_thermalized(v_f):
    return get_norm(v_f) < THERMALIZED_VELOCITY_THRESHOLD


def simulate_simple(x, v):
    while True:
        distance = get_distance_to_next_interaction(MACROSCOPIC_CROSS_SECTION)
        x = x + ((distance / get_norm(v)) * v)

        if is_outside_right(x):
            return "ESCAPED_RIGHT", x
        if is_outside_left(x):
            return "ESCAPED_LEFT", x
        if is_absorbed():
            return "ABSORBED", x

        v, theta_lab, energy = elastic_collision(v, get_atomic_mass_target())

        if is_thermalized(v):
            return "THERMALIZED", x



