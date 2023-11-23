from math import floor


lst = [-35, -11, -4, 0, 6, 21, 43]

old_min = min(lst)
old_range = max(lst) - old_min

new_min = 0
new_range = 255 - new_min

norm = [(n - old_min) / old_range * new_range + new_min for n in lst]
print(f"{lst}\n{norm}")
