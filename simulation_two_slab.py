import math
from random import random

from numpy.linalg import norm as get_norm

from main import elastic_collision

MACROSCOPIC_CS_ABSORPTION_WATER = 0.010063
MACROSCOPIC_CS_ABSORPTION_CARBON = 0.00026
MACROSCOPIC_CS_SCATTERING_WATER = 0.361663 - 0.010063
MACROSCOPIC_CS_SCATTERING_CARBON = 0.3846

HYDROGEN_SCATTERING_PERCENTAGE = 7.8 / 10.5
ABSORBANCE_RATIO_WATER = MACROSCOPIC_CS_ABSORPTION_WATER / 0.361663
ABSORBANCE_RATIO_CARBON = MACROSCOPIC_CS_ABSORPTION_CARBON / (0.3846 + MACROSCOPIC_CS_ABSORPTION_CARBON)
# Unitless, initial velocity 1MeV, thermalized velocity threshold 1eV
THERMALIZED_VELOCITY_THRESHOLD = 1 / math.sqrt(10 ** 6)
# In cm-1
MACROSCOPIC_CROSS_SECTION_WATER = 0.361663
MACROSCOPIC_CROSS_SECTION_CARBON = 0.3846 + MACROSCOPIC_CS_ABSORPTION_CARBON  # cm-1


def get_media(x, water_width):
    position = x[0]
    if position < 0 or position > 30:
        return "VOID"
    elif position < water_width:
        return "WATER"
    else:
        return "CARBON"


def get_macroscopic_cross_section(current_media):
    if current_media == "WATER":
        return MACROSCOPIC_CROSS_SECTION_WATER
    else:
        # Using the same cross-section for carbon and void, but it is not relevant.
        return MACROSCOPIC_CROSS_SECTION_CARBON


def get_macroscopic_cross_section_absorbance(current_media):
    if current_media == "WATER":
        return MACROSCOPIC_CS_ABSORPTION_WATER
    else:
        # Using the same cross-section for carbon and void, but it is not relevant.
        return MACROSCOPIC_CS_ABSORPTION_CARBON


def get_macroscopic_cross_section_scattering(current_media):
    if current_media == "WATER":
        return MACROSCOPIC_CS_SCATTERING_WATER
    else:
        # Using the same cross-section for carbon and void, but it is not relevant.
        return MACROSCOPIC_CS_SCATTERING_CARBON


def get_other_macroscopic_cross_section(current_media):
    if current_media == "WATER":
        return MACROSCOPIC_CROSS_SECTION_CARBON
    else:
        return MACROSCOPIC_CROSS_SECTION_WATER


def get_next_position(current_position, v, water_width):
    current_media = get_media(current_position, water_width)
    macroscopic_cross_section = get_macroscopic_cross_section(current_media)
    minus_log_r = - math.log(random())
    distance_to_next_interaction = minus_log_r / macroscopic_cross_section
    next_position = current_position + ((distance_to_next_interaction / get_norm(v)) * v)

    next_media = get_media(next_position, water_width)

    if next_media == current_media:
        return next_position
    if v[0] > 0 and current_position[0] > water_width:  # Escaped right
        return next_position
    if v[0] < 0 and current_position[0] < water_width:  # Escaped left
        return next_position
    else:  # Media change from water to graphite or the other way around.
        next_macroscopic_cs = get_other_macroscopic_cross_section(current_media)

        position_change_medium = current_position + v * ((water_width - current_position[0]) / v[0])
        # print("####################### Start #####################")
        # print(f"Current position is : {current_position}")
        # print(f"Previous next position is : {next_position}")
        # print(f"The position_change_medium is : {position_change_medium}")
        # print(f"Velocity is : {v}")

        d1 = get_norm(position_change_medium - current_position)
        # print(f"d1 : {d1}")
        d2 = (minus_log_r - d1 * macroscopic_cross_section) / next_macroscopic_cs
        # print(f"d2 : {d2}")
        dt = d1 + d2
        # print(f"dt : {dt}")
        # print(f"Previous dt : {distance_to_next_interaction}")

        next_position = current_position + ((dt / get_norm(v)) * v)
        # print(f"Recalculated next position is : {next_position}")

        return next_position


def get_atomic_mass_target(x, water_width):
    media = get_media(x, water_width)

    if media == "WATER":
        if random() < HYDROGEN_SCATTERING_PERCENTAGE:
            return 1
        else:
            return 16
    elif media == "CARBON":
        return 12


def is_outside_right(position):
    return position[0] > 30


def is_outside_left(position):
    return position[0] < 0


def is_absorbed(x, water_width):
    media = get_media(x, water_width)
    if media == "WATER":
        return random() < ABSORBANCE_RATIO_WATER
    elif media == "CARBON":
        return random() < ABSORBANCE_RATIO_CARBON


def is_thermalized(v_f):
    return get_norm(v_f) < THERMALIZED_VELOCITY_THRESHOLD


def simulate_simple_two_slab(x, v, water_width):
    while True:
        x = get_next_position(x, v, water_width)

        if is_outside_right(x):
            return "ESCAPED_RIGHT", x
        if is_outside_left(x):
            return "ESCAPED_LEFT", x
        if is_absorbed(x, water_width):
            return "ABSORBED", x

        v, theta_lab, energy = elastic_collision(v, get_atomic_mass_target(x, water_width))

        if is_thermalized(v):
            return "THERMALIZED", x


def simulate_single_collision_two_slab(x, v, water_width):
    x = get_next_position(x, v, water_width)

    if is_outside_right(x):
        return "ESCAPED_RIGHT", x
    if is_absorbed(x, water_width):
        return "ABSORBED", x
    return "OTHER", x


def simulate_multiple_collision_two_slab(x, v, water_width):
    collisions = 0
    while True:
        x = get_next_position(x, v, water_width)

        if is_outside_right(x):
            return "ESCAPED_RIGHT", x, collisions
        if is_outside_left(x):
            return "ESCAPED_LEFT", x, collisions

        collisions += 1

        if is_absorbed(x, water_width):
            return "ABSORBED", x, collisions

        v, theta_lab, energy = elastic_collision(v, get_atomic_mass_target(x, water_width))

        if is_thermalized(v):
            return "THERMALIZED", x, collisions
