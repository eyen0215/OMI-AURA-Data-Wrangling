import h5py
import numpy as np
import matplotlib.pyplot as plt
import imageio

def to_gif(arrays, name, cmap = 'viridis'):
    # Create a list to store image frames for the GIF
    frames = []

    # Plot each array and convert it to an image
    for index, array in enumerate(arrays):
        # Plot the array using a colormap
        fig, ax = plt.subplots()
        ax.imshow(array, cmap=cmap)  
        ax.axis('off')

        # Add the frame number as text on the plot
        ax.text(5, 5, f'Frame {index + 1}', color='white', fontsize=12, 
                bbox=dict(facecolor='black', alpha=0.7))

        # Convert the plot to a NumPy array
        fig.canvas.draw()
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

        # add the image to the list
        frames.append(image)

        plt.close(fig)

    # Create a GIF using the collected frames with looping enabled and slower animation
    output_gif = f'src/{name}'
    imageio.mimsave(output_gif, frames, fps=2, loop=0)  # Set fps to a lower value for a slower animation

    print(f"GIF saved as {output_gif}")


# Open the HDF5 file and read the arrays
sulfur_file = 'Sulfur_arrays.h5'
loss_file = 'Loss_arrays.h5'

sulfur_arrays = []
loss_arrays = []

# Read all datasets from the HDF5 files
with h5py.File(sulfur_file, 'r') as f, h5py.File(loss_file, 'r') as g:
    for key in f.keys():
        sulfur_arrays.append(f[key][:])  
    for key in g.keys():
        loss_arrays.append(g[key][:])

# Combine the sulfur and loss arrays for plotting
to_gif(sulfur_arrays, 'Sulfur_array_animation.gif', 'viridis')
to_gif(loss_arrays, 'Loss_array_animation.gif', 'viridis')


