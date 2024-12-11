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


def get_atomic_mass_target(medium):
    if medium == "OXYGEN":
        return 16
    elif medium == "HYDROGEN":
        return 1
    elif random() < HYDROGEN_SCATTERING_PERCENTAGE:
        return 1
    else:
        return 16


def is_absorbed():
    return random() < ABSORBANCE_PERCENTAGE


def is_thermalized(v_f):
    return get_norm(v_f) < THERMALIZED_VELOCITY_THRESHOLD


def simulate_water_simple_with_absorbance(v):
    num_collisions = 0
    medium = "WATER"
    while True:
        num_collisions += 1
        if is_absorbed():
            return "ABSORBED", num_collisions

        v, theta_lab, energy = elastic_collision(v, get_atomic_mass_target(medium))

        if is_thermalized(v):
            return "THERMALIZED", num_collisions


def simulate_homogenous_without_absorbance(v, medium):
    num_collisions = 0
    while True:
        num_collisions += 1
        v, theta_lab, energy = elastic_collision(v, get_atomic_mass_target(medium))

        if is_thermalized(v):
            return num_collisions
