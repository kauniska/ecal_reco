'''
Geometric properties of the detector
'''

n_strips = 24
n_layers = 8

#Dimensions of the cells expressed in cm
"""
# without the correct geometry implemented
width = 1
thickness = 1
thickness_screen = 0
"""
width = 1.6
thickness = 0.4
thickness_screen = 0.6

#clockcycle value
clockcycle_value = 6.25 #nanosecond

#Height of the detector
total_height = 2*n_layers*(thickness+thickness_screen)

#Distance between 2 active layer of the same orientation
Delta_z = 2*(thickness+thickness_screen)

# #General offset on time computed with time alignement algorithm, for each channel
# muX = np.nan_to_num(np.ndarray(shape=(8,64), dtype=float), nan=0, posinf=0, neginf=0)*0
# muY = np.nan_to_num(np.ndarray(shape=(8,64), dtype=float), nan=0, posinf=0, neginf=0)*0


## Delay of each SiPM channel (in picosecond)
SiPM_delay = [870.358, 901.456, 902.071, 265.938, 768.12, 1022.583, 1021.051, 392.153, 653.118, 1145.177, 1143.645, 513.848, 551.309, 1266.305, 1264.772, 634.464, 436.307, 1388.898, 1387.366, 756.116, 334.069, 1509.722, 1508.493, 876.691, 786.982, 770.555, 766.403, 110.931, 665.725, 890.342, 886.62, 229.49, 545.865, 1012.454, 1008.303, 351.174, 424.537, 1132.241, 1128.519, 471.39, 304.677, 1254.218, 1250.202, 593.073, 129.977, 1372.852, 1370.418, 713.289, 715.466, 837.501, 838.851, 715.666, 594.212, 956.807, 958.157, 593.983, 473.567, 1077.149, 1078.929, 473.338, 351.883, 1196.455, 1198.234, 352.513, 220.503, 1317.091, 1319.006, 232.297, 103.387, 1430.947, 1438.312, 113.068, 966.758, 940.508, 944.237, 312.913, 841.689, 1059.079, 1081.92, 433.354, 703.427, 1177.281, 1202.692, 555.037, 578.358, 1368.033, 1321.998, 675.253, 440.095, 1488.376, 1442.77, 796.251, 296.295, 1602.636, 1562.075, 915.877]