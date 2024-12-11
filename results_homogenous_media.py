import math
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from simulation import simulate_simple, simulate_single_collision, simulate_multiple_collision
from simulation_homogenous_media import simulate_homogenous_without_absorbance, simulate_water_simple_with_absorbance

MACROSCOPIC_CS_ABSORBANCE = 0.010063


def print_average_collision_number(num_simulations, media):
    collisions_count = defaultdict(int)

    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        collisions = simulate_homogenous_without_absorbance(initial_velocity, media)

        collisions_count[collisions] += 1

    total_collisions = sum(key * value for key, value in collisions_count.items())
    total_simulations = sum(collisions_count.values())

    average_collisions = total_collisions / total_simulations if total_simulations > 0 else 0

    print(f"Average number of collisions in {media}: {average_collisions}")


def print_average_collision_number_with_absorbance(num_simulations):
    collisions_count = defaultdict(int)

    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, collisions = simulate_water_simple_with_absorbance(initial_velocity)
        if result == "THERMALIZED":
            collisions_count[collisions] += 1

    total_collisions = sum(key * value for key, value in collisions_count.items())
    total_simulations = sum(collisions_count.values())

    average_collisions = total_collisions / total_simulations if total_simulations > 0 else 0

    print(f"Average number of collisions in water with absorbance: {average_collisions}")

simulation_number=100000
print_average_collision_number(simulation_number,"HYDROGEN")
print_average_collision_number(simulation_number,"OXYGEN")
print_average_collision_number(simulation_number,"WATER")
print_average_collision_number_with_absorbance(simulation_number)
