import numpy as np 
import matplotlib.pyplot as plt
from simulation import simulate

ENERGY_ANALYSED_COLLISIONS = 8

# Initial conditions
initial_position = np.array([0, 0, 0])
initial_velocity = np.array([1, 0, 0])

# Simulate neutron transport and collect results
num_simulations = 1000
termination_counts = {"ESCAPED_RIGHT": 0, "ESCAPED_LEFT": 0, "ABSORBED": 0, "THERMALIZED": 0}
thermalized_positions = []
scattering_angles = []
collision_energies = [[] for i in range(ENERGY_ANALYSED_COLLISIONS)]

# Average Energy Theory
ratio_H = 1/2
ratio_O = 257/289
ratio_H2O = (2*3.9*ratio_H + 1*2.7*ratio_O) / (2*3.9 + 1*2.7)
n_collisions = int(-np.log(10**6) / np.log(ratio_H2O))
Ps = 0.3516 / (0.3516 + 0.010063)
Ptherm = Ps**n_collisions
alpha_H = 0
alpha_O = (15/17)**2
alpha_H2O = (2*3.9*alpha_H + 1*2.7*alpha_O) / (2*3.9 + 1*2.7)

cos_H = 1/np.sqrt(2)
cos_O = 1/np.sqrt(1+16**2)
cos_H2O = (2*3.9*cos_H + 1*2.7*cos_O) / (2*3.9 + 1*2.7)

# Store the energies of the first five neutrons for the last plot
first_five_neutron_energies = []

for i in range(num_simulations):
    result, positions, angles, energies = simulate(initial_position, initial_velocity, ENERGY_ANALYSED_COLLISIONS=8)
    termination_counts[result] += 1
    
    if result == "THERMALIZED":
        thermalized_positions.append(positions[-1][0])  # Collect x-coordinate only
    
    scattering_angles.extend(angles)
    
    for j in range(min(len(energies), ENERGY_ANALYSED_COLLISIONS)):
        collision_energies[j].append(energies[j])
    
    # Store energies of the first five neutrons for plotting
    if i < 5:
        first_five_neutron_energies.append(energies)


# Calculer la position moyenne de thermalisation
mean_thermalization_position = np.mean(thermalized_positions)

# Afficher le résultat
print("Position moyenne de thermalisation :", mean_thermalization_position)





# Plot histogram of thermalized positions
plt.hist(thermalized_positions, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Position where neutrons are thermalized (cm)')
plt.ylabel('Count')
plt.title('Distribution of Thermalized Positions')
plt.show()



# Plot histogram of scattering angles with added mean and cos_H2O lines
plt.hist(scattering_angles, bins=50, color='salmon', edgecolor='black')
plt.xlabel('Scattering angle (radians)')
plt.ylabel('Count')
plt.title('Distribution of Scattering Angles in LAB')

# Calculate the mean scattering angle
mean_scattering_angle = np.mean(scattering_angles)

# Add a thick solid green line for cos_H2O angle
plt.axvline(np.arccos(cos_H2O), color='green', linestyle='-', linewidth=4, label='Theoretical Mean Angle')

# Add a blue dotted line for the mean angle
plt.axvline(mean_scattering_angle, color='blue', linestyle=':', linewidth=2, label='Mean Angle')

plt.legend()
plt.show()


# Plot histogram of scattering angles
plt.hist(np.cos(scattering_angles), bins=50, color='salmon', edgecolor='black')
plt.xlabel('cos (scattering angle)')
plt.ylabel('Count')
plt.title('Cosine of the Scattering Angle in LAB for Hydrogen')
plt.show()



# Plotting the energy evolution for the first five neutrons
plt.figure()
for i, energies in enumerate(first_five_neutron_energies):
    plt.plot(range(1, len(energies) + 1), energies, marker='o', label=f'Neutron {i+1}')
    
plt.xlabel("Collision Number")
plt.ylabel("Energy (MeV)")
plt.title("Evolution of Energy over Eight First Collisions for Five Neutrons")
plt.legend()
plt.show()


# Plot histogram for energy frequency after each of the first 8 collisions
for i in range(ENERGY_ANALYSED_COLLISIONS):
    plt.hist(collision_energies[i], bins=50, color='lightgreen', edgecolor='black')
    plt.xlabel('Energy after collision {} (MeV)'.format(i+1))
    plt.ylabel('Count')
    plt.title(f'Energy Distribution After Collision {i+1}')
    
    if i == 0:
        plt.axvline(alpha_O, color='black', linestyle='--', linewidth=1.5, label='E = αO')
        plt.legend()
    
    plt.show()
    

plt.hist(collision_energies[0], bins=50, color='lightgreen', edgecolor='black')
plt.xlabel('Energy after collision {} (MeV)'.format(1))
plt.ylabel('Count')
plt.title('Energy Distribution After Collision 1')

plt.show()    


# Calculate the average energy per collision over all simulations and plot it
average_energy_per_collision = [
    np.mean(collision_energies[i]) for i in range(ENERGY_ANALYSED_COLLISIONS)
]

# Calculate theoretical energy values based on the number of collisions
theoretical_energy_per_collision = [
    1 * (ratio_H2O ** n) for n in range(1, ENERGY_ANALYSED_COLLISIONS + 1)
]

plt.plot(range(1, ENERGY_ANALYSED_COLLISIONS + 1), average_energy_per_collision, marker='o', linestyle='-', color='blue', label='Simulated Average Energy')
plt.plot(range(1, ENERGY_ANALYSED_COLLISIONS + 1), theoretical_energy_per_collision, marker='x', linestyle='--', color='red', label='Theoretical Energy')
plt.xlabel('Collision Number')
plt.ylabel('Energy (MeV)')
plt.title('Evolution of Average Energy over Eight First Collisions')
plt.legend()
plt.show()

plt.plot(range(1, ENERGY_ANALYSED_COLLISIONS + 1), average_energy_per_collision, marker='o', linestyle='-', color='blue')
plt.xlabel('Collision Number')
plt.ylabel('Average Energy (MeV)')
plt.title('Evolution of Average Energy over Eight first Collisions')
plt.show()



# Plot the proportion of each termination type
plt.bar(termination_counts.keys(), termination_counts.values(), color='lightcoral', edgecolor='black')
plt.xlabel('Termination Type')
plt.ylabel('Count')
plt.title('Proportion of Each Termination Type')
plt.show()


# Calculate normalized proportions of absorbed and thermalized neutrons
total_absorbed_thermalized = termination_counts["ABSORBED"] + termination_counts["THERMALIZED"]
proportion_absorbed = termination_counts["ABSORBED"] / total_absorbed_thermalized
proportion_thermalized = termination_counts["THERMALIZED"] / total_absorbed_thermalized

plt.bar([' Absorbed', ' Thermalized'], 
        [proportion_absorbed, proportion_thermalized], 
        color='blue', label='Simulated', alpha=0.6)
plt.bar(['Absorbed', 'Thermalized'], 
        [1 - Ptherm, Ptherm], 
        color='orange', label='Theoretical', alpha=0.6)
plt.ylabel('Proportion')
plt.title('Comparison of Simulated vs Theoretical Proportions of Absorbed and Thermalized Neutrons')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))  # Moves legend outside the plot
plt.tight_layout()  # Adjusts plot to fit everything nicely
plt.show()