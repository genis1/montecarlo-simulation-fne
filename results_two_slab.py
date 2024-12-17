import math
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from simulation_two_slab import simulate_simple_two_slab, get_media, \
    get_macroscopic_cross_section_absorbance, simulate_single_collision_two_slab, simulate_multiple_collision_two_slab, \
    get_macroscopic_cross_section_scattering


def plot_flux_over_position(num_simulations, num_buckets, water_width, color):
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
        cs_absorbance = get_macroscopic_cross_section_absorbance(
            get_media([position + (bucket_width / 2), 0, 0], water_width))
        return count * initial_neutron_flux / (
                num_simulations * bucket_width * cs_absorbance)

    flux_dict = {(bucket * bucket_width + bucket_width / 2): transform(bucket * bucket_width, count) for bucket, count
                 in
                 bucket_counts.items()}
    flux_dict = {key: value for key, value in sorted(flux_dict.items())}

    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.plot(positions, flux, color=color, linestyle='-', linewidth=2, label=f"Water {water_width}cm")
    plt.xlabel('Distance (cm)')
    plt.ylabel('Flux ($cm^{-2}s^{-1}$)')
    plt.title('Flux over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    plt.axvline(x=water_width, color=color, linestyle='--', linewidth=0.5)
    print("Simulation done")


def plot_flux_over_position_times_4(num_simulations_1, num_buckets_1):
    plot_flux_over_position(num_simulations=num_simulations_1, num_buckets=num_buckets_1, water_width=30,
                            color='skyblue')
    plot_flux_over_position(num_simulations=num_simulations_1, num_buckets=num_buckets_1, water_width=5, color='brown')
    plot_flux_over_position(num_simulations=num_simulations_1, num_buckets=num_buckets_1, water_width=10,
                            color='olive')
    plot_flux_over_position(num_simulations=num_simulations_1, num_buckets=num_buckets_1, water_width=15,
                            color='limegreen')
    plt.legend()
    file_name = f"total_flux_{num_simulations_1}_num_buckets{num_buckets_1}"
    plt.savefig(f"figures/{file_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()


def plot_single_collision_flux(num_simulations, num_buckets, water_width, color):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position = simulate_single_collision_two_slab(initial_position, initial_velocity, water_width)
        if result != "OTHER":
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    def transform(position, count):
        cs_scattering = get_macroscopic_cross_section_scattering(
            get_media([position + (bucket_width / 2), 0, 0], water_width))
        return count * initial_neutron_flux / (
                num_simulations * bucket_width * cs_scattering)

    flux_dict = {(bucket * bucket_width + (bucket_width / 2)): transform(bucket * bucket_width, count) for bucket, count
                 in
                 bucket_counts.items()}
    flux_dict = {key: value for key, value in sorted(flux_dict.items())}

    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.plot(positions, flux, color=color, linestyle='-', linewidth=0.5, label=f"Water {water_width}cm")
    plt.xlabel('Distance (cm)')
    plt.ylabel('Flux ($cm^{-2}s^{-1}$)')
    plt.title('Flux over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    plt.axvline(x=water_width, color=color, linestyle='--', linewidth=0.5)
    print("Simulation done")


def plot_single_collision_flux_times_4(num_simulations, num_buckets):
    plot_single_collision_flux(num_simulations=num_simulations, num_buckets=num_buckets, water_width=30,
                               color='skyblue')
    plot_single_collision_flux(num_simulations=num_simulations, num_buckets=num_buckets, water_width=5, color='brown')
    plot_single_collision_flux(num_simulations=num_simulations, num_buckets=num_buckets, water_width=10,
                               color='olive')
    plot_single_collision_flux(num_simulations=num_simulations, num_buckets=num_buckets, water_width=15,
                               color='limegreen')
    plt.legend()
    file_name = f"single_collision_flux_{num_simulations}_num_buckets{num_buckets}"
    plt.savefig(f"figures/{file_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()


def plot_multiple_collision_flux_over_position(num_simulations, num_buckets, water_width, color):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position, collisions = simulate_multiple_collision_two_slab(initial_position, initial_velocity,
                                                                                  water_width)
        if result != "ABSORBED" or collisions == 1:
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    def transform(position, count):
        cs_absorbance = get_macroscopic_cross_section_absorbance(
            get_media([position + (bucket_width / 2), 0, 0], water_width))
        return count * initial_neutron_flux / (
                num_simulations * bucket_width * cs_absorbance)

    flux_dict = {(bucket * bucket_width + (bucket_width / 2)): transform(bucket * bucket_width, count) for bucket, count
                 in
                 bucket_counts.items()}
    flux_dict = {key: value for key, value in sorted(flux_dict.items())}

    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.plot(positions, flux, color=color, linestyle='-', linewidth=2, label=f"Water {water_width}cm")
    plt.xlabel('Distance (cm)')
    plt.ylabel('Flux ($cm^{-2}s^{-1}$)')
    plt.title('Flux over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    plt.axvline(x=water_width, color=color, linestyle='--', linewidth=0.5)
    print("Simulation done")


def plot_multiple_collision_flux_over_position_times_4(num_simulations, num_buckets):
    plot_multiple_collision_flux_over_position(num_simulations=num_simulations, num_buckets=num_buckets, water_width=30,
                                               color='skyblue')
    plot_multiple_collision_flux_over_position(num_simulations=num_simulations, num_buckets=num_buckets, water_width=5,
                                               color='brown')
    plot_multiple_collision_flux_over_position(num_simulations=num_simulations, num_buckets=num_buckets, water_width=10,
                                               color='olive')
    plot_multiple_collision_flux_over_position(num_simulations=num_simulations, num_buckets=num_buckets, water_width=15,
                                               color='limegreen')
    plt.legend()
    file_name = f"multiple_collision_flux_{num_simulations}_num_buckets{num_buckets}"
    plt.savefig(f"figures/{file_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()


def plot_slowing_down_density(num_simulations, num_buckets, water_width, color):
    bucket_width = 30 / num_buckets  # cm
    bucket_counts = defaultdict(int)
    initial_neutron_flux = 1000  # n/cm^2 s

    initial_position = np.array([0, 0, 0])
    initial_velocity = np.array([1, 0, 0])

    for i in range(num_simulations):
        result, final_position = simulate_simple_two_slab(initial_position, initial_velocity, water_width)
        if result != "THERMALIZED":
            continue

        final_x = final_position[0]

        bucket_number = math.floor(final_x / bucket_width)
        bucket_counts[bucket_number] += 1

    def transform(count):
        return count * initial_neutron_flux / (
                num_simulations * bucket_width)

    flux_dict = {(bucket * bucket_width + (bucket_width / 2)): transform(count) for bucket, count in
                 bucket_counts.items()}
    flux_dict = {key: value for key, value in sorted(flux_dict.items())}

    flux = np.array(list(flux_dict.values()))
    positions = np.array(list(flux_dict.keys()))

    plt.plot(positions, flux, color=color, linestyle='-', linewidth=2, label=f"Water {water_width}cm")
    plt.xlabel('Distance (cm)')
    plt.ylabel('Slowing down density ($cm^{-3}s^{-1}$)')
    plt.title('Slowing down density over position for an initial flux of 1000$cm^{-2}s^{-1}$')
    plt.axvline(x=water_width, color=color, linestyle='--', linewidth=0.5)
    print("Simulation done")



def plot_slowing_down_density_times_4(num_simulations, num_buckets):
    plot_slowing_down_density(num_simulations=num_simulations, num_buckets=num_buckets, water_width=30,
                              color='skyblue')
    plot_slowing_down_density(num_simulations=num_simulations, num_buckets=num_buckets, water_width=5, color='brown')
    plot_slowing_down_density(num_simulations=num_simulations, num_buckets=num_buckets, water_width=10,
                              color='olive')
    plot_slowing_down_density(num_simulations=num_simulations, num_buckets=num_buckets, water_width=15,
                              color='limegreen')
    plt.legend()
    file_name = f"slowing_down_density_{num_simulations}_num_buckets{num_buckets}"
    plt.savefig(f"figures/{file_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()


scale_num_simulation = 1  # 100 for standard, 1000 for good results, 1 for quick test
os.makedirs("figures", exist_ok=True)

plot_flux_over_position_times_4(2000 * scale_num_simulation, 60)
print("Figure 1")

plot_single_collision_flux_times_4(num_simulations=50000 * scale_num_simulation, num_buckets=600)
print("Figure 2")

plot_multiple_collision_flux_over_position_times_4(num_simulations=2000 * scale_num_simulation, num_buckets=60)
print("Figure 3")

plot_slowing_down_density_times_4(num_simulations=2000 * scale_num_simulation, num_buckets=60)
print("Figure 4")
