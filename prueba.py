import numpy as np

# Cucl4:
# [1.86062,
#  0.00219, 0.03274, 0.00015, 0.00006,
#  0.00219, 0.03274, 0.00015, 0.00006,
#  0.00219, 0.03274, 0.00015, 0.00006,
#  0.00219, 0.03274, 0.00015, 0.00006,
#  ]

# Cu(mnt)2:


suma = 0
for i in range(0, len(density)):
    suma += density[i]

metal = float(density[0])
metal_density = metal * 100 / suma
print(np.round(metal_density, 1))