import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap


def read_design_file(filename):
    """Read a design file and return x positions and heights"""
    data = np.loadtxt(filename, comments='#')
    if data.size == 0:
        return np.array([]), np.array([])

    if data.ndim == 1:  # Only one data point
        x = np.array([data[0]])
        h = np.array([data[1]])
    else:
        x = data[:, 0]
        h = data[:, 1]

    return x, h


def plot_all_iterations():
    """Plot all iterations on the same figure"""
    files = sorted(glob.glob('design_iter_*.txt'))

    if not files:
        print("No design files found!")
        return

    plt.figure(figsize=(12, 8))

    # Create a colormap for iterations
    colors = plt.cm.viridis(np.linspace(0, 1, len(files)))

    for i, file in enumerate(files):
        x, h = read_design_file(file)
        if len(x) > 0:
            iteration = int(file.split('_')[-1].split('.')[0])
            plt.plot(x,
                     h,
                     color=colors[i],
                     label=f'Iter {iteration}',
                     alpha=0.7,
                     linewidth=2)

    plt.xlabel('Position (x)', fontsize=14)
    plt.ylabel('Beam Height (h)', fontsize=14)
    plt.title('Beam Design Evolution', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('plots/design_evolution_all.png', dpi=300, bbox_inches='tight')


def plot_iteration_comparison():
    """Plot specific iterations for comparison"""
    files = sorted(glob.glob('design_iter_*.txt'))

    if not files:
        print("No design files found!")
        return

    # Select specific iterations to compare
    iterations_to_plot = [0, 10, 50, 100, -1]  # First, some middle, and last

    plt.figure(figsize=(12, 8))

    colors = ['red', 'blue', 'green', 'orange', 'purple']
    line_styles = ['-', '--', '-.', ':', '-']

    for i, iter_idx in enumerate(iterations_to_plot):
        if iter_idx == -1:
            file = files[-1]  # Last iteration
        else:
            # Find file with this iteration number
            matching_files = [
                f for f in files if f.endswith(f'{iter_idx:06d}.txt')
            ]
            if not matching_files:
                continue
            file = matching_files[0]

        x, h = read_design_file(file)
        if len(x) > 0:
            iteration = int(file.split('_')[-1].split('.')[0])
            plt.plot(x,
                     h,
                     color=colors[i % len(colors)],
                     linestyle=line_styles[i % len(line_styles)],
                     label=f'Iter {iteration}',
                     linewidth=3,
                     alpha=0.8)

    plt.xlabel('Position (x)', fontsize=14)
    plt.ylabel('Beam Height (h)', fontsize=14)
    plt.title('Beam Design - Key Iterations', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/design_key_iterations.png', dpi=300)


def create_animation():
    """Create an animation of the design evolution"""
    files = sorted(glob.glob('design_iter_*.txt'))

    if not files:
        print("No design files found!")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    # Read first file to get x values
    x, _ = read_design_file(files[0])

    line, = ax.plot(x, np.zeros_like(x), 'b-', linewidth=2)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, 1)  # Adjust based on your design variable range
    ax.set_xlabel('Position (x)', fontsize=12)
    ax.set_ylabel('Beam Height (h)', fontsize=12)
    ax.set_title('Beam Design Evolution', fontsize=14)
    ax.grid(True, alpha=0.3)

    def update(frame):
        file = files[frame]
        x, h = read_design_file(file)
        iteration = int(file.split('_')[-1].split('.')[0])
        line.set_ydata(h)
        ax.set_title(f'Beam Design - Iteration {iteration}', fontsize=14)
        return line,

    ani = FuncAnimation(fig,
                        update,
                        frames=len(files),
                        interval=200,
                        blit=True)

    # Save animation
    ani.save('design_evolution.gif', writer='pillow', fps=5)
    print("Animation saved as design_evolution.gif")

    plt.close()


def plot_convergence():
    """Plot the maximum and minimum heights over iterations"""
    files = sorted(glob.glob('design_iter_*.txt'))

    if not files:
        print("No design files found!")
        return

    iterations = []
    max_heights = []
    min_heights = []
    mean_heights = []

    for file in files:
        x, h = read_design_file(file)
        if len(h) > 0:
            iteration = int(file.split('_')[-1].split('.')[0])
            iterations.append(iteration)
            max_heights.append(h.max())
            min_heights.append(h.min())
            mean_heights.append(h.mean())

    plt.figure(figsize=(12, 8))

    plt.plot(iterations, max_heights, 'r-', label='Max Height', linewidth=2)
    plt.plot(iterations, min_heights, 'b-', label='Min Height', linewidth=2)
    plt.plot(iterations, mean_heights, 'g--', label='Mean Height', linewidth=2)

    plt.xlabel('Iteration', fontsize=14)
    plt.ylabel('Height', fontsize=14)
    plt.title('Beam Height Statistics Over Iterations', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/height_statistics.png', dpi=300)


if __name__ == "__main__":

    if len(os.sys.argv) > 1:
        choice = os.sys.argv[1]
    else:
        print("Beam Design Visualization")
        print("1. Plot all iterations")
        print("2. Plot key iterations comparison")
        print("3. Create animation")
        print("4. Plot height statistics")
        choice = input("Enter your choice (1-4): ").strip()

    if os.path.exists('plots') is False:
        os.makedirs('plots')

    if choice == '1':
        plot_all_iterations()
    elif choice == '2':
        plot_iteration_comparison()
    elif choice == '3':
        create_animation()
    elif choice == '4':
        plot_convergence()
    elif choice.lower() == 'all':
        print("Running all visualizations...")
        plot_all_iterations()
        plot_iteration_comparison()
        create_animation()
        plot_convergence()
