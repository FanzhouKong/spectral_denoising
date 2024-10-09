import matplotlib.pyplot as plt
import numpy as np
import spectral_denoising as sd
# x = np.linspace(0, 2 * np.pi, 200)
# y = np.sin(x)

# fig, ax = plt.subplots()
# ax.plot(x, y)
# plt.show()
peak = np.array([[48.992496490478516 ,154.0],
                  [63.006099700927734, 265.0],
                  [79.02062225341797, 521.0],
                  [159.02373146795, 999]], 
                  
                  dtype = np.float32)

sd.head_to_tail_plot(peak, peak, pmz = 200)
plt.show()
