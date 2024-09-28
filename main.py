import random
import math


def elastic_collision_with_stationary_mass(v0, m0, m2):
    # theta is the angle of deflection with respect to the original trajectory, 0 for no deflection
    cos_theta = 2 * random.random() - 1
    # phi being the angle of direction of the deflection
    phi = 2 * math.pi * random.random()
    # Calculating v0
    alpha = (m2 ** 2 + 2 * m2 * m0 * cos_theta + m0 ** 2) / ((m2 + m0) ** 2)
    print("alpha: " + str(alpha))
    v0f = v0 * pow(alpha, 0.5)
    # angle definition in equations:
    sin_theta = pow(1 - cos_theta ** 2, 0.5)
    vx = v0f * cos_theta
    vy = v0f * sin_theta * math.sin(phi)
    vz = v0f * sin_theta * math.cos(phi)
    return [vx, vy, vz]


def sum_of_squares(arr):
    return sum(x * x for x in arr)


v0 = [100, 0, 0]
m0 = 1
m2 = 100

print("v0: " + str(v0))

E0 = m0 * pow(v0[0], 2) / 2
print("E0: " + str(E0))
vf = elastic_collision_with_stationary_mass(v0[0], m0, m2)
print("vf: " + str(vf))

Ef = sum_of_squares(vf) * m0 / 2
print("Ef: " + str(Ef))
