import os

file_path = 'resource/LaPalma1971_slope_profile.txt'

with open(file_path, 'r') as f:
    lines = f.readlines()

# Get the last line
last_line = lines[-1].strip().split('\t')
last_dist = float(last_line[3])

# Calculate average slope of the last 5 points
slopes = [float(l.strip().split('\t')[4]) for l in lines[-5:]]
avg_slope = sum(slopes) / len(slopes)

# Extend up to 3000m
new_lines = []
current_dist = last_dist + 10.0
while current_dist <= 3000.0:
    new_lines.append(f"0\t0\t0\t{int(current_dist)}\t{avg_slope}\n")
    current_dist += 10.0

with open(file_path, 'a') as f:
    for line in new_lines:
        f.write(line)

print(f"Extended {file_path} to 3000m with average slope {avg_slope:.2f}")
