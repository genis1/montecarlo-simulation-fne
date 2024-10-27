import numpy as np
import matplotlib.pyplot as plt
from simulation import simulate

ENERGY_ANALYSED_COLLISIONS = 8

# Initial conditions
initial_position = np.array([0, 0, 0])
initial_velocity = np.array([1, 0, 0])

# Simulate neutron transport and collect results
num_simulations = 10000
termination_counts = {"ESCAPED_RIGHT": 0, "ESCAPED_LEFT": 0, "ABSORBED": 0, "THERMALIZED": 0}
thermalized_positions = []
scattering_angles = []
collision_energies = [[] for i in range(ENERGY_ANALYSED_COLLISIONS)]

for i in range(num_simulations):
    result, positions, angles, energies = simulate(initial_position, initial_velocity, ENERGY_ANALYSED_COLLISIONS = 8)
    termination_counts[result] += 1
    
    if result == "THERMALIZED":
        thermalized_positions.append(positions[-1][0])  # Collect x-coordinate only
        
    scattering_angles.extend(angles)
    
    for j in range(min(len(energies), ENERGY_ANALYSED_COLLISIONS)):
        collision_energies[j].append(energies[j])

# Plot histogram of thermalized positions
plt.hist(thermalized_positions, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Position where neutrons are thermalized (cm)')
plt.ylabel('Frequency')
plt.title('Distribution of Thermalized Positions')
plt.show()

# Plot histogram of scattering angles
plt.hist(scattering_angles, bins=50, color='salmon', edgecolor='black')
plt.xlabel('Scattering angle (radians)')
plt.ylabel('Frequency')
plt.title('Distribution of Scattering Angles in LAB')
plt.show()

# Plot histogram for energy frequency after each of the first 8 collisions
for i in range(ENERGY_ANALYSED_COLLISIONS):
    plt.hist(collision_energies[i], bins=50, color='lightgreen', edgecolor='black')
    plt.xlabel('Energy after collision {} (MeV)'.format(i+1))
    plt.ylabel('Frequency')
    plt.title(f'Energy Distribution After Collision {i+1}')
    plt.show()
    
# Plot the proportion of each termination type
plt.bar(termination_counts.keys(), termination_counts.values(), color='lightcoral', edgecolor='black')
plt.xlabel('Termination Type')
plt.ylabel('Count')
plt.title('Proportion of Each Termination Type')
plt.show()

# Calculate the average energy per collision over all simulations and plot it
average_energy_per_collision = [
    np.mean(collision_energies[i]) for i in range(ENERGY_ANALYSED_COLLISIONS)
]

plt.plot(range(1, ENERGY_ANALYSED_COLLISIONS + 1), average_energy_per_collision, marker='o', linestyle='-', color='blue')
plt.xlabel('Collision Number')
plt.ylabel('Average Energy (eV)')
plt.title('Evolution of Average Energy over Eight first Collisions')
plt.show()
