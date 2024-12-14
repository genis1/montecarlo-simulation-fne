import math
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from simulation import simulate_simple, simulate_single_collision, simulate_multiple_collision
from simulation_two_slab import MACROSCOPIC_CS_SCATTERING_WATER

MACROSCOPIC_CS_ABSORBANCE = 0.010063


def plot_flux_over_position(num_simulations, num_buckets):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position = simulate_simple(initial_position, initial_velocity)
        if result != "ABSORBED":
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    def transform(count):
        return count * initial_neutron_flux / (
                num_simulations * bucket_width * MACROSCOPIC_CS_ABSORBANCE)

    flux_dict = {(bucket * bucket_width + bucket_width / 2): transform(count) for bucket, count
                 in
                 bucket_counts.items()}
    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.bar(positions, flux, color='skyblue', edgecolor='black', width=bucket_width)
    plt.xlabel('Distance (cm)')
    plt.ylabel('Flux ($cm^{-2}s^{-1}$)')
    plt.title('Flux over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    file_name = f"total_flux_{num_simulations}_num_buckets{num_buckets}"
    plt.savefig(f"figures/task1/{file_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()


def plot_single_collision_flux(num_simulations, num_buckets):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position = simulate_single_collision(initial_position, initial_velocity)
        if result != "OTHER":
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    def transform(count):
        return count * initial_neutron_flux / (
                num_simulations * bucket_width * MACROSCOPIC_CS_SCATTERING_WATER)

    flux_dict = {(bucket * bucket_width + (bucket_width / 2)): transform(count) for bucket, count
                 in
                 bucket_counts.items()}
    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.bar(positions, flux, color='skyblue', edgecolor='skyblue', width=bucket_width)
    plt.xlabel('Distance (cm)')
    plt.ylabel('Flux ($cm^{-2}s^{-1}$)')
    plt.title('Singe collision flux over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    file_name = f"single_collision_flux_{num_simulations}_num_buckets{num_buckets}"
    plt.savefig(f"figures/task1/{file_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()


def plot_multiple_collision_flux_over_position(num_simulations, num_buckets):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position, collisions = simulate_multiple_collision(initial_position, initial_velocity)
        if result != "ABSORBED" or collisions == 1:
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    def transform(count):
        return count * initial_neutron_flux / (
                num_simulations * bucket_width * MACROSCOPIC_CS_ABSORBANCE)

    flux_dict = {(bucket * bucket_width + (bucket_width / 2)): transform(count) for bucket, count
                 in
                 bucket_counts.items()}
    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.bar(positions, flux, color='skyblue', edgecolor='black', width=bucket_width)
    plt.xlabel('Distance (cm)')
    plt.ylabel('Flux ($cm^{-2}s^{-1}$)')
    plt.title('Multiple collision flux over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    file_name = f"multiple_collision_flux_{num_simulations}_num_buckets{num_buckets}"
    plt.savefig(f"figures/task1/{file_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()


def plot_slowing_down_density(num_simulations, num_buckets):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position = simulate_simple(initial_position, initial_velocity)
        if result != "THERMALIZED":
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    def transform(count):
        return count * initial_neutron_flux / (num_simulations * bucket_width)

    flux_dict = {(bucket * bucket_width + (bucket_width / 2)): transform(count) for bucket, count in
                 bucket_counts.items()}

    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.bar(positions, flux, color='skyblue', edgecolor='black', width=bucket_width)
    plt.xlabel('Distance (cm)')
    plt.ylabel('Slowing down density ($cm^{-3}s^{-1}$)')
    plt.title('Slowing down density over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    file_name = f"slowing_down_density_{num_simulations}_num_buckets{num_buckets}"
    plt.savefig(f"figures/task1/{file_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()


os.makedirs("figures", exist_ok=True)
os.makedirs("figures/task1", exist_ok=True)

scale_num_simulation = 1000  # 100 for standard, 1000 for good results, 1 for quick test
plot_flux_over_position(num_simulations=200 * scale_num_simulation, num_buckets=100)
plot_single_collision_flux(num_simulations=5000 * scale_num_simulation, num_buckets=600)
plot_multiple_collision_flux_over_position(num_simulations=200 * scale_num_simulation, num_buckets=100)
plot_slowing_down_density(num_simulations=200 * scale_num_simulation, num_buckets=60)
