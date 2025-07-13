def mhz_to_hz(mhz):
    return mhz * 1e+6

def cycles_to_seconds(cycles, frequency):
    return cycles / frequency

k = 5
max_mhz = 4100
min_mhz = 800
max_hz = mhz_to_hz(max_mhz)
min_hz = mhz_to_hz(min_mhz)
num_pixels = 5000
max_iter = 50
cycles_per_distance = 50
cycles_per_average = 100
minutes = 90
sampled_fps = 2

cycles_average_per_frame = cycles_per_average * k * max_iter
cycles_distance_per_frame = cycles_per_distance * num_pixels * k * max_iter

frames = minutes * 60 * sampled_fps

total_cycles = (cycles_average_per_frame + cycles_distance_per_frame) * frames
total_time_max = cycles_to_seconds(total_cycles, min_hz)
total_time_min = cycles_to_seconds(total_cycles, max_hz)

print("Assuming KMeans is the slowest thread, a", minutes, "minute video will",
      "take between", total_time_min / 60, "and", total_time_max / 60, "minutes.")
