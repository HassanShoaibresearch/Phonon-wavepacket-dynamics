import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage import uniform_filter1d

# Directory containing the data files
data_folder = "data3"

# Get all file names in the data folder
file_names = os.listdir(data_folder)

# Filter and sort the file names based on time
file_names = sorted([file_name for file_name in file_names if file_name.endswith(".dump")], key=lambda x: float(x.split(".")[0]))

# Initialize lists to store time and energy data
time_data = []
energy_data = []

# Read data from each file
for file_name in file_names:
    # Extract the time value from the file name
    time = float(file_name.split(".")[0])
    time_data.append(time)

    # Read the data from the file
    file_path = os.path.join(data_folder, file_name)
    k_data, energy = np.loadtxt(file_path, unpack=True)

    energy_data.append(energy)

# Convert energy_data to a numpy array
energy_data = np.array(energy_data)

# Calculate the average energy for each k
average_energies = np.mean(energy_data[9000:], axis=0)

# Calculate the total energy for each k
total_energy = np.sum(energy_data, axis=1)

# Find the indices of the top 10 average energy values
top_indices = np.argsort(average_energies)[-10:]

# Initialize lists to store the selected k and energy values
selected_ks = []
selected_energies = []

# Get the selected k and energy values
for index in top_indices:
    selected_ks.append(k_data[index])
    selected_energies.append(average_energies[index])

# Assuming energy_data is your original array of size (40, 200)
# Assuming selected_ks is your list of size 10 containing the indices of selected K values

# Convert selected_ks to a NumPy array and ensure it has integer data type
selected_ks = np.array(selected_ks, dtype=np.int)

# Use NumPy indexing to select the desired columns from energy_data
selected_energy_data = energy_data[:, selected_ks]

# The resulting array, selected_energy_data, will have size (40, 10)

# Plot the data
#E_smooth= gaussian_filter1d(selected_energy_data, sigma=2)
# Smoothing the selected_energy_data using a moving average
E_smooth = uniform_filter1d(selected_energy_data, size=15, axis=0, mode='nearest')

# Plotting the selected_energy_data
for i in range(selected_energy_data.shape[1]):
    label = "K{}".format(selected_ks[i])  # Label based on selected_ks array
    plt.plot(time_data, selected_energy_data[:, i], label=label)

# Add labels and legend
plt.xlabel('Time')
plt.ylabel('Energy')
plt.legend()
plt.grid()

# Display the plot
plt.show()

# Plotting the selected_energy_data
for i in range(selected_energy_data.shape[1]):
    label = "K{}".format(selected_ks[i])  # Label based on selected_ks array
    plt.plot(time_data, selected_energy_data[:, i], label=label)

# Add labels and legend
plt.xlabel('Time')
plt.ylabel('Energy (log)')
plt.legend()
plt.yscale('log')
plt.grid()

# Display the plot
plt.show()
    
    # Plotting the selected_energy_data
for i in range(selected_energy_data.shape[1]):
    label = "K{}".format(selected_ks[i])  # Label based on selected_ks array
    plt.plot(time_data, E_smooth[:, i], label=label)

# Add labels and legend
plt.xlabel('Time')
plt.ylabel('Energy (log-smoothed)')
plt.legend()
plt.yscale('log')
plt.grid()
# Display the plot
plt.show()

# Plotting the total_energy
# Smoothing the selected_energy_data using a moving average
E_total_smooth = uniform_filter1d(total_energy, size=10, axis=0, mode='nearest')
plt.plot(time_data, E_total_smooth)

# Add labels and legend
plt.xlabel('Time')
plt.ylabel('Total Energy (smoothed)')
plt.legend()
plt.yscale('log')
plt.grid()
# Display the plot
plt.show()