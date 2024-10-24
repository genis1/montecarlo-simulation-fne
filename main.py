import math
import random

import numpy as np
from numpy.linalg import norm as get_norm


def elastic_collision(v_0_LAB, mass_target):
    cos_theta_collision_CM = 2 * random.random() - 1

    Acos = mass_target * cos_theta_collision_CM
    common_numerator = mass_target + 1
    to_be_squared = math.pow(mass_target, 2) - math.pow(Acos, 2)

    parallel_factor = (Acos + 1) / common_numerator
    perpendicular_factor = get_norm(v_0_LAB) * math.sqrt(to_be_squared) / common_numerator

    v_f_perpendicular_LAB = normalize_vector(np.cross(v_0_LAB, get_random_vector()))

    v_f_LAB = parallel_factor * v_0_LAB + perpendicular_factor * v_f_perpendicular_LAB

    return v_f_LAB


def normalize_vector(vector):
    norm = get_norm(vector)
    if norm == 0:
        return vector
    return vector / norm


def get_random_vector():
    random_cos_theta = 2 * random.random() - 1
    random_sin_theta = math.sqrt(1 - math.pow(random_cos_theta, 2))
    random_phi = 2 * math.pi * random.random()
    random_cos_phi = math.cos(random_phi)
    random_sin_phi = math.sin(random_phi)
    random_vector = [
        random_cos_phi * random_sin_theta,
        random_sin_phi * random_sin_theta,
        random_cos_theta
    ]
    return np.array(random_vector)


v2 = np.array([1 / 2, 1 / 2, 1 / 2])

v_f = elastic_collision(v2, 0.1)
print("elastic collsion: " + str(v_f))
print("norm vf:" + str(get_norm(v_f)))

min_values = min(get_norm(elastic_collision(v2, 16)) for x in range(1000))
print("min_value" + str(min_values))
