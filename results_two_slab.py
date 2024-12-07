import math
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from simulation_two_slab import simulate_simple_two_slab, get_media, \
    get_macroscopic_cross_section_absorbance, simulate_single_collision_two_slab


def plot_flux_over_position(num_simulations, num_buckets, water_width):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position = simulate_simple_two_slab(initial_position, initial_velocity, water_width)
        if result != "ABSORBED":
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    def transform(position, count):
        cs_absorbance = get_macroscopic_cross_section_absorbance(get_media([position, 0, 0], water_width))
        return count * initial_neutron_flux / (
                num_simulations * bucket_width * cs_absorbance)

    flux_dict = {bucket * bucket_width: transform(bucket * bucket_width, count) for bucket, count in
                 bucket_counts.items()}
    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.bar(positions, flux, color='skyblue', edgecolor='black', width=bucket_width)
    plt.xlabel('Distance (cm)')
    plt.ylabel('Flux ($cm^{-2}s^{-1}$)')
    plt.title('Flux over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    plt.axvline(x=water_width, color='red', linestyle='--', linewidth=1.5)
    plt.show()

def plot_single_collision_flux(num_simulations, num_buckets, water_width):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position = simulate_single_collision_two_slab(initial_position, initial_velocity, water_width)
        if result != "ABSORBED":
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    def transform(position, count):
        cs_absorbance = get_macroscopic_cross_section_absorbance(get_media([position, 0, 0], water_width))
        return count * initial_neutron_flux / (
                num_simulations * bucket_width * cs_absorbance)

    flux_dict = {bucket * bucket_width: transform(bucket * bucket_width, count) for bucket, count in
                 bucket_counts.items()}
    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.bar(positions, flux, color='skyblue', edgecolor='skyblue', width=bucket_width)
    plt.xlabel('Distance (cm)')
    plt.ylabel('Flux ($cm^{-2}s^{-1}$)')
    plt.title('Flux over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    plt.axvline(x=water_width, color='red', linestyle='--', linewidth=1.5)
    plt.show()


plot_flux_over_position(num_simulations=200000, num_buckets=60, water_width=5)
plot_flux_over_position(num_simulations=200000, num_buckets=60, water_width=10)
plot_flux_over_position(num_simulations=200000, num_buckets=60, water_width=15)
plot_single_collision_flux(num_simulations=5000000, num_buckets=500,water_width=5)
plot_single_collision_flux(num_simulations=5000000, num_buckets=500,water_width=10)
plot_single_collision_flux(num_simulations=5000000, num_buckets=500,water_width=15)
