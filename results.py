import math
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from simulation import simulate_absorption

MACROSCOPIC_CS_ABSORBANCE = 0.010063


def plot_flux_over_position(num_simulations, num_buckets):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position = simulate_absorption(initial_position, initial_velocity)
        if result != "ABSORBED":
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    buckets = np.array(list(bucket_counts.keys())) * bucket_width
    counts = np.array(list(bucket_counts.values())) * initial_neutron_flux / (
            num_simulations * bucket_width * MACROSCOPIC_CS_ABSORBANCE)
    plt.bar(buckets, counts, color='skyblue', edgecolor='black', width=bucket_width)
    plt.xlabel('Distance (cm)')
    plt.ylabel('Flux ($cm^{-2}s^{-1}$)')
    plt.title('Flux over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    plt.show()


plot_flux_over_position(num_simulations=100000, num_buckets=100)