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


def simulate(x, v, ENERGY_ANALYSED_COLLISIONS):
    positions = []
    scattering_angles = []
    energies = []

    num_collisions = 0

    while True:
        distance = get_distance_to_next_interaction(MACROSCOPIC_CROSS_SECTION)
        x = x + ((distance / get_norm(v)) * v)

        if is_outside_right(x):
            return "ESCAPED_RIGHT", positions, scattering_angles, energies, num_collisions
        if is_outside_left(x):
            return "ESCAPED_LEFT", positions, scattering_angles, energies, num_collisions
        if is_absorbed():
            return "ABSORBED", positions, scattering_angles, energies, num_collisions

        v, theta_lab, energy = elastic_collision(v, get_atomic_mass_target())
        scattering_angles.append(theta_lab)

        if num_collisions < ENERGY_ANALYSED_COLLISIONS:
            energies.append(energy)
        num_collisions += 1
        if is_thermalized(v):
            positions.append(x)
            return "THERMALIZED", positions, scattering_angles, energies, num_collisions

        positions.append(x)


def simulate_absorption(x, v):
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


def simulate_single_collision(x, v):
    distance = get_distance_to_next_interaction(MACROSCOPIC_CROSS_SECTION)
    x = x + ((distance / get_norm(v)) * v)

    if is_outside_right(x):
        return "ESCAPED_RIGHT", x
    if is_absorbed():
        return "ABSORBED", x
    return "OTHER", x
