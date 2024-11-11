import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def initialize_lattice(size):
    # Create a 2D array with dimensions size x size
    lattice = np.zeros((size, size), dtype=int)
    return lattice


def compute_neighbor_indices(size):
    # Create an empty dictionary to hold neighbor indices
    neighbor_indices = {}
    
    # Iterate over each position in the lattice
    for x in range(size):
        for y in range(size):
            # Calculate neighbor positions with wrapping
            neighbors = [((x - 1) % size, y),  
                        ((x + 1) % size, y),  
                        (x, (y - 1) % size),  
                        (x, (y + 1) % size)]
            # Store neighbors in the dictionary
            neighbor_indices[(x, y)] = neighbors
    return neighbor_indices
    

def calculate_interaction_energy(lattice, site, particle, neighbor_indices, epsilon_AA, epsilon_BB, epsilon_AB):
    # Get coordinates of the site
    x, y = site
    interaction_energy = 0
    
    # Iterate over each neighbor
    for neighbor in neighbor_indices[(x, y)]:
        neighbor_particle = lattice[neighbor[0]][neighbor[1]]
        
        # Check if the neighbor is occupied
        if neighbor_particle != 0:
            if particle == 1: 
                if neighbor_particle == 1: 
                    interaction_energy += epsilon_AA
                else: 
                    interaction_energy += epsilon_AB
            else:  
                if neighbor_particle == 2:  
                    interaction_energy += epsilon_BB
                else: 
                    interaction_energy += epsilon_AB
    return interaction_energy


def attempt_move(lattice, N_A, N_B, N_empty, neighbor_indices, params):
    size = lattice.shape[0]  # Dimension of lattice
    N_sites = size * size
    beta = 1 / params['T']
    
    # Extract parameters
    epsilon_A = params['epsilon_A']
    epsilon_B = params['epsilon_B']
    epsilon_AA = params['epsilon_AA']
    epsilon_BB = params['epsilon_BB']
    epsilon_AB = params['epsilon_AB']
    mu_A = params['mu_A']
    mu_B = params['mu_B']
    
    # Decide whether to add or remove a particle (50% chance each)
    if np.random.rand() < 0.5: #Adding
        if N_empty == 0:
            return N_A, N_B, N_empty 
        
        # Select a random site from empty sites
        empty_sites = np.argwhere(lattice == 0)
        site = empty_sites[np.random.choice(empty_sites.shape[0])]
        
        # Decide which particle to add (A or B) with equal probability
        if np.random.rand() < 0.5:  
            particle = 1
            mu = mu_A
            epsilon = epsilon_A
            N_s = N_A
        else: 
            particle = 2
            mu = mu_B
            epsilon = epsilon_B
            N_s = N_B
        
        # Calculate delta_E
        delta_E = epsilon + calculate_interaction_energy(lattice, tuple(site), particle, neighbor_indices, epsilon_AA, epsilon_BB, epsilon_AB)
        
        # Calculate acceptance probability
        acc_prob = min(1, (N_empty) / (N_s + 1) * np.exp(-beta * (delta_E - mu)))
        
        # Generate a random number
        r = np.random.rand()
        if r < acc_prob:

            lattice[site[0], site[1]] = particle
            if particle == 1:
                N_A += 1
            else:
                N_B += 1
            N_empty -= 1
    else:  # Removing
        if N_sites - N_empty == 0:
            return N_A, N_B, N_empty  # No particles to remove
        
        # Select a random site from occupied sites
        occupied_sites = np.argwhere(lattice != 0)
        site = occupied_sites[np.random.choice(occupied_sites.shape[0])]
        
        # Get the particle at site
        particle = lattice[site[0], site[1]]
        
        if particle == 1:
            mu = mu_A
            epsilon = epsilon_A
            N_s = N_A
        else: 
            mu = mu_B
            epsilon = epsilon_B
            N_s = N_B
        
        # Calculate delta_E
        delta_E = -epsilon - calculate_interaction_energy(lattice, tuple(site), particle, neighbor_indices, epsilon_AA, epsilon_BB, epsilon_AB)
        
        # Calculate acceptance probability
        acc_prob = min(1, N_s / (N_empty + 1) * np.exp(-beta * (delta_E + mu)))
        
        # Generate a random number
        r = np.random.rand()
        if r < acc_prob:
            lattice[site[0], site[1]] = 0  # Remove particle
            if particle == 1:
                N_A -= 1
            else:
                N_B -= 1
            N_empty += 1
            
    return N_A, N_B, N_empty


def run_simulation(size, n_steps, params):
    # Initialize the lattice and compute neighbor indices
    lattice = initialize_lattice(size)
    neighbor_indices = compute_neighbor_indices(size)
    N_sites = size * size
    
    # Initialize counts
    N_A = 0
    N_B = 0
    N_empty = N_sites
    
    # Create arrays for coverage
    coverage_A = np.zeros(n_steps)
    coverage_B = np.zeros(n_steps)

    for step in range(n_steps):
        N_A, N_B, N_empty = attempt_move(lattice, N_A, N_B, N_empty, neighbor_indices, params)
        
        # Update coverage
        coverage_A[step] = N_A / N_sites
        coverage_B[step] = N_B / N_sites

    return lattice, coverage_A, coverage_B


def plot_lattice(lattice, ax, title):
    size = lattice.shape[0]  # Get the dimension of the lattice


    for x in range(size):
        for y in range(size):
            if lattice[x, y] == 1:  # Particle A
                ax.plot(x + 0.5, y + 0.5, 'ro', markersize=10)  
            elif lattice[x, y] == 2:  # Particle B
                ax.plot(x + 0.5, y + 0.5, 'bo', markersize=10)  

    # Set axis limits and labels
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_xticks(np.arange(0, size + 1, 1))  # Keep the ticks
    ax.set_yticks(np.arange(0, size + 1, 1))  # Keep the ticks
    ax.set_xticklabels(['']*(size+1))  # Remove x tick labels
    ax.set_yticklabels(['']*(size+1))  # Remove y tick labels
    ax.grid()
    ax.set_title(title, fontsize=10)  # Set the title of the axis

    return ax

#Example Simulation to Test Implementation
size = 4
n_steps = 10000
mus_A = np.linspace(-0.2, 0, 7)
Ts = np.linspace(0.001, 0.019, 7)
params = []
for mu_A in mus_A:
    for T in Ts:
        params.append({'epsilon_A': -0.1,
        'epsilon_B': -0.1,
        'epsilon_AA': 0,
        'epsilon_BB': 0,
        'epsilon_AB': 0,
        'mu_A': mu_A,
        'mu_B': -0.1,
        'T': T})
# Run the simulation
np.random.seed(42)
final_lattice = np.zeros((len(mus_A), len(Ts), size, size))
mean_coverage_A = np.zeros((len(mus_A), len(Ts)))
mean_coverage_B = np.zeros((len(mus_A), len(Ts)))
for i, param in enumerate(params):
    lattice, coverage_A, coverage_B = run_simulation(size, n_steps, param)
    final_lattice[i // len(Ts), i % len(Ts)] = lattice
    mean_coverage_A[i // len(Ts), i % len(Ts)] = np.mean(coverage_A[-1000:])
    mean_coverage_B[i // len(Ts), i % len(Ts)] = np.mean(coverage_B[-1000:])

# Plot the T-mu_A phase diagram
fig, axs = plt.subplot_mosaic([[0, 1, 2], [3, 4, 5]], figsize=(6.5, 4.5))

# Mean coverage of A
axs[0].pcolormesh(mus_A, Ts, mean_coverage_A.T, cmap='viridis', vmin=0, vmax=1)
axs[0].set_title(r'$\langle \theta_A \rangle$')
axs[0].set_xlabel(r'$\mu_A$ (eV)')
axs[0].set_ylabel(r'T (K)')
axs[0].set_yticks([Ts[0], Ts[1], Ts[2], Ts[3], Ts[4], Ts[5], Ts[6]])
axs[0].set_yticklabels([round(Ts[0]/0.0000861733), round(Ts[1]/0.0000861733), round(Ts[2]/0.0000861733), round(Ts[3]/0.0000861733), round(Ts[4]/0.0000861733), round(Ts[5]/0.0000861733), round(Ts[6]/0.0000861733)])

# Mean coverage of B
axs[1].pcolormesh(mus_A, Ts, mean_coverage_B.T, cmap='viridis', vmin=0, vmax=1)
axs[1].set_title(r'$\langle \theta_B \rangle$')
axs[1].set_xlabel(r'$\mu_A$ (eV)')
axs[1].set_yticks([])

# Mean total coverage
cax = axs[2].pcolormesh(mus_A, Ts, mean_coverage_A.T + mean_coverage_B.T, cmap='viridis', vmin=0, vmax=1)
axs[2].set_title(r'$\langle \theta_A + \theta_B \rangle$')
axs[2].set_xlabel(r'$\mu_A$ (eV)')
axs[2].set_yticks([])
fig.colorbar(cax, ax=axs[2], location='right', fraction=0.1)

# Plot the final lattice configuration

# mu_A = -0.2 eV and T = 116K
axs[3] = plot_lattice(final_lattice[0, 3], axs[3], r'$\mu_A = -0.2$ eV, $T = 116K$')

# mu_A = -0.1 eV and T = 116K
axs[4] = plot_lattice(final_lattice[3, 3], axs[4], r'$\mu_A = -0.1$ eV, $T = 116K$')

# mu_A = 0 eV and T = 116K
axs[5] = plot_lattice(final_lattice[6, 3], axs[5], r'$\mu_A = 0$ eV, $T = 116K$')

plt.tight_layout()
plt.show()
plt.savefig("Example_Simulation.png")