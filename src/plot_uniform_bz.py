import numpy as np
import h5py
import matplotlib.pyplot as plt

# Load the HDF5 file
with h5py.File("uniform_bz.h5", "r") as f:
    print("Available keys:", list(f.keys()))

    # Read datasets
    x  = f["x"][:]
    vx = f["vx"][:]
    vy = f["vy"][:]
    vz = f["vz"][:]

# Plot velocity components as a function of position
plt.figure(figsize=(10, 6))
plt.plot(x, vx, label="vx")
plt.plot(x, vy, label="vy")
plt.plot(x, vz, label="vz")

plt.xlabel("x (position)")
plt.ylabel("velocity components")
plt.title("Velocity vs Position from uniform_bz.h5")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

