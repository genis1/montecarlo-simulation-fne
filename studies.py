import numpy as np
import matplotlib.pyplot as plt
from simulation import simulate

# Initial conditions
initial_position = np.array([0, 0, 0])
initial_velocity = np.array([1, 0, 0])

# Simulate neutron transport and collect results
num_simulations = 1000
results = []
thermalized_positions = []
scattering_angles = []

for i in range(num_simulations):
    result, positions, angles = simulate(initial_position, initial_velocity)
    results.append(result)
    
    if result == "THERMALIZED":
        thermalized_positions.append(positions[-1][0])  # Collect x-coordinate only
    scattering_angles.extend(angles)

# Calculate average scattering angle
average_angle = np.mean(scattering_angles)

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

# Print average scattering angle
print(f'Average scattering angle in LAB: {average_angle:.2f} radians')