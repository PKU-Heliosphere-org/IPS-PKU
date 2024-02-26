#
# geometer.py
# Geometer routines.
# Coded by Dr. L. Zhang <lifeaint1987@hotmail.com>.


maxloc = lambda x: x.index(max(x))
circle = lambda x: x[1:] + [x[0]]
rewind_at_i = lambda l, i: l[i + 1:] + l[0:i + 1]


def rewind_point_list(list_x, list_y):
    if (len(list_x) < 3 or len(list_y) < 3):
        return (list_x, list_y)
    distance = lambda l: (l[0] - l[1]) ** 2 + (l[2] - l[3]) ** 2
    distances = list(map(distance, zip(list_x, circle(list_x), list_y, circle(list_y))))
    i_max = maxloc(distances)
    if (i_max == len(distances) - 1):
        return (list_x, list_y)
    else:
        return (rewind_at_i(list_x, i_max), rewind_at_i(list_y, i_max))
