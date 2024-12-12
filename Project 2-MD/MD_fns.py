import matplotlib.pyplot as plt
import numpy as np

#Simulation Parameters
k_B = 1.38e-23 #J/K
dt = 0.01  # Time step
total_steps = 10000  # Number of steps
box_size = 100.0  # Size of the cubic box
k = 1.0  # Spring constant
mass = 1.0  # Particle mass
r0 = 1.0  # Equilibrium bond length
target_temperature = 0.1  # Target temperature
rescale_interval = 100  # Steps between velocity rescaling
n_particles = 20  # Number of particles
epsilon_repulsive = 1.0  # Depth of repulsive LJ potential
epsilon_attractive = 0.5  # Depth of attractive LJ potential
sigma = 1.0  # LJ potential parameter
cutoff = 2**(1/6) * sigma
position_trajectory = []
velocity_trajectory = []
forces_trajectory = []
temperature_trajectory = []



#Functions
def apply_pbc(position, box_size):
    #Applies periodic boundary conditions
    return position % box_size


def initialize_chain(n_particles, box_size, r0):
    positions = np.zeros((n_particles, 3))

    # Initialize the first position at the center of the box
    current_position = np.array([box_size/2, box_size/2, box_size/2])
    positions[0] = current_position

    # Iterate to initialize the rest of the positions in a chain
    for i in range(1, n_particles):
        direction = np.random.randn(3) 
        direction /= np.linalg.norm(direction)  # Normalize the vector to unit length
        next_position = current_position + r0 * direction #Move direction from current position
        positions[i] = apply_pbc(next_position, box_size)
        current_position = positions[i] #Update current position
    return positions

def initialize_velocities(n_particles, target_temperature, mass):
    velocities = np.random.randn(n_particles, 3) 
    velocities *= np.sqrt(k_B * target_temperature / mass) #Random velocities put to Maxwell-Boltzmann Dist
    velocities -= np.mean(velocities, axis=0) #Remove net momentum
    return velocities

def minimum_image(displacement, box_size):
    #Applies minimum image to account for PBC
    for i in range(3): #Adjusts displacement if outside box
        if displacement[i] > box_size / 2: 
            displacement[i] -= box_size
        elif displacement[i] < -box_size / 2:
            displacement[i] += box_size
    return displacement    

def compute_harmonic_forces(positions, k, r0, box_size):
    forces = np.zeros_like(positions) #Intialize force array
    for i in range(0, (n_particles-2)):
        displacement = positions[i+1] - positions[i]
        displacement = minimum_image(displacement, box_size) #Apply PBC to displacement vector
        distance = np.linalg.norm(displacement)
        force_magnitude = -k * (distance - r0) #Harmonic force magnitude
        force = force_magnitude * (displacement / distance)

        #apply forces to both molecules in opposite directions
        forces[i] -= force
        forces[i+1] += force
    return forces

def compute_lennard_jones_forces(positions, epsilon, sigma, box_size, interaction_type):
    forces = np.zeros_like(positions)
    for i in range(0, (n_particles-1)):
        for j in range((i+1), (n_particles-1)): #Determines interaction type based on given conditions
            if interaction_type == 'repulsive' and np.linalg.norm(i-j) == 2:
                epsilon = float(epsilon_repulsive) #given
            elif interaction_type == 'attractive' and np.linalg.norm(i-j) > 2:
                epsilon = float(epsilon_attractive) #given
            else:
                continue
            displacement = positions[j] - positions[i]
            displacement = minimum_image(displacement, box_size) #Apply PBC to displacement vector
            distance = float(np.linalg.norm(displacement))

            if distance < cutoff:
                lj = (sigma / distance)**12 - (.5 * (sigma / distance)**6)
                epsilon = float(epsilon)
                force_magnitude = 24 * epsilon * (lj) / distance
                force = force_magnitude * (displacement / distance) #Normalize force vector

                #Apply forces to both molecules in opposite directions
                forces[i] -= force
                forces[j] += force 
    return forces

def velocity_verlet(positions, velocities, forces, dt, mass):
    #Update velocities and positions
    velocities += 0.5 * forces / mass * dt
    positions += velocities * dt
    positions = apply_pbc(positions, box_size) #boundary conditions

    #New forces and velocities from updated positions
    forces_harmonic = compute_harmonic_forces(positions, k, r0, box_size)
    forces_repulsive = compute_lennard_jones_forces(positions, epsilon_repulsive, sigma, box_size, 'repulsive')
    forces_attractive = compute_lennard_jones_forces(positions, epsilon_attractive, sigma, box_size, 'attractive')
    forces_new = forces_harmonic + forces_repulsive + forces_attractive
    velocities += 0.5 * forces_new / mass * dt
    return positions, velocities, forces_new

def rescale_velocities(velocities, target_temperature, mass):
    kinetic_energy = 0.5 * mass * np.sum(np.linalg.norm(velocities, axis=1)**2)
    current_temperature = (2 / 3) * kinetic_energy / (n_particles * k_B) #calculate current temp from KE
    scaling_factor = np.sqrt(target_temperature / current_temperature)
    velocities *= scaling_factor #scales currect velocity based on target and current temps
    return velocities




#Example Simulation

# Initialize positions and velocities
positions = initialize_chain(n_particles, box_size, r0)
velocities = initialize_velocities(n_particles, target_temperature, mass)

# Simulation loop
for step in range(total_steps):
    # Compute forces
    forces_harmonic = compute_harmonic_forces(positions, k, r0, box_size)
    forces_repulsive = compute_lennard_jones_forces(positions, epsilon_repulsive, sigma, box_size, 'repulsive')
    forces_attractive = compute_lennard_jones_forces(positions, epsilon_attractive, sigma, box_size, 'attractive')
    total_forces = forces_harmonic + forces_repulsive + forces_attractive
    
    # Integrate equations of motion
    positions, velocities, total_forces = velocity_verlet(positions, velocities, total_forces, dt, mass)
    
    # Apply thermostat
    if step % rescale_interval == 0:
        velocities = rescale_velocities(velocities, target_temperature, mass)

# Store data for analysis
    if step % 100 == 0:  # Store every 100 steps
        position_trajectory.append(positions.copy())  # Store a copy of the positions
        velocity_trajectory.append(velocities.copy())  # Store a copy of the velocities
        forces_trajectory.append(total_forces.copy())  # Store a copy of total forces


#Analysis

def calculate_radius_of_gyration(positions):
    center_of_mass = np.mean(positions, axis = 0)
    Rg_squared = np.mean(sum((positions - center_of_mass)^2, axis = 1))
    Rg = np.sqrt(Rg_squared)
    
    return Rg

def calculate_end_to_end_distance(positions):
    Ree = np.linalg.norm(positions[-1] - posiions[0])
    
    return Ree


#Example Analysis

# Arrays to store properties
temperatures = np.linspace(0.1, 1.0, 10)
Rg_values = []
Ree_values = []
potential_energies = []

for T in temperatures:
    # Set target temperature
    target_temperature = T
    # (Re-initialize positions and velocities)
    # (Run simulation)
    # Compute properties
    Rg = calculate_radius_of_gyration(positions)
    Ree = calculate_end_to_end_distance(positions)
    Rg_values.append(Rg)
    Ree_values.append(Ree)
    potential_energies.append(np.mean(potential_energy_array))

# Plotting
plt.figure()
plt.plot(temperatures, Rg_values, label='Radius of Gyration')
plt.xlabel('Temperature')
plt.ylabel('Radius of Gyration')
plt.title('Radius of Gyration vs Temperature')
plt.legend()
plt.show()

plt.figure()
plt.plot(temperatures, Ree_values, label='End-to-End Distance')
plt.xlabel('Temperature')
plt.ylabel('End-to-End Distance')
plt.title('End-to-End Distance vs Temperature')
plt.legend()
plt.show()

plt.figure()
plt.plot(temperatures, potential_energies, label='Potential Energy')
plt.xlabel('Temperature')
plt.ylabel('Potential Energy')
plt.title('Potential Energy vs Temperature')
plt.legend()
plt.show()