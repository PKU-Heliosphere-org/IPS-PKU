#
# geometer.py
# Geometer routines.
# Coded by Dr. L. Zhang <lifeaint1987@hotmail.com>.

import numpy as np


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


cross = lambda a, b: [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
]

def local_frame_unit_vector(theta0, phi0):
    e_zeta = [np.sin(theta0) * np.cos(phi0), np.sin(theta0) * np.sin(phi0), np.cos(theta0)]
    e_ksi = [np.sin(theta0 + np.pi/2) * np.cos(phi0), np.sin(theta0 + np.pi/2) * np.sin(phi0), np.cos(theta0 + np.pi/2)]
    e_eta = cross(e_zeta, e_ksi)
    return (e_ksi, e_eta, e_zeta)


def local_ksi_eta_zeta_to_xyz(vector_ksi, theta0, phi0):
    e_ksi, e_eta, e_zeta = local_frame_unit_vector(theta0, phi0)
    prod = lambda i: vector_ksi[0] * e_ksi[i] + vector_ksi[1] * e_eta[i] + vector_ksi[2] * e_zeta[i]
    return (prod(0), prod(1), prod(2))


def global_theta_phi_from_xyz(loc):
    r = (loc[0] ** 2 + loc[1] ** 2 + loc[2] ** 2) ** 0.5
    theta = np.arccos(loc[2] / r)
    phi = np.arctan2(loc[1], loc[0]) # Numpy has arctan2(y, x) to compute arg(x + yi).
    if (phi < 0):
        phi += 2 * np.pi
    return (theta, phi)


def az_alt_deg_to_theta_phi_rad(az_deg, alt_deg):
    theta_deg = 90 - alt_deg
    phi_deg = 360 - az_deg
    theta_rad, phi_rad = np.deg2rad(theta_deg), np.deg2rad(phi_deg)
    return (theta_rad, phi_rad)


def theta_phi_rad_to_az_alt_deg(theta_rad, phi_rad):
    theta_deg, phi_deg = np.rad2deg(theta_rad), np.rad2deg(phi_rad)
    az_deg = 360 - phi_deg
    alt_deg = 90 - theta_deg
    return (az_deg, alt_deg)


def compute_one_point_az_alt(local_theta_rad, local_phi_rad, theta0_rad, phi0_rad):
    loc_ksi = [
        np.sin(local_theta_rad) * np.cos(local_phi_rad),
        np.sin(local_theta_rad) * np.sin(local_phi_rad),
        np.cos(local_theta_rad)
    ]
    loc_xyz = local_ksi_eta_zeta_to_xyz(loc_ksi, theta0_rad, phi0_rad)
    global_theta_rad, global_phi_rad = global_theta_phi_from_xyz(loc_xyz)
    return theta_phi_rad_to_az_alt_deg(global_theta_rad, global_phi_rad)


def compute_cme_circle(az_alt_deg_sun_viewed, r_cme):
    theta0, phi0 = az_alt_deg_to_theta_phi_rad(az_alt_deg_sun_viewed[0], az_alt_deg_sun_viewed[1])
    local_theta_rad = 2 * np.arcsin(r_cme * np.pi / 180 * 1000 / 3600 / 2)
    local_phi_rad = np.arange(start=0, stop=np.pi * 2 + np.pi * 2 / 72, step=np.pi * 2 / 72)
    return [compute_one_point_az_alt(local_theta_rad, the_phi, theta0, phi0) for the_phi in local_phi_rad]


def example_self_test():
    print(compute_cme_circle([30, 60], 10))
    return True


# Cut the x_list, y_list in TWO parts if possible to cut.
def cut_x_y_lists(x_list, y_list):
    threshold = 0.25 # as r * r. Hence the difference will be 0.5.
    distance = lambda l: (l[0] - l[1]) ** 2 + (l[2] - l[3]) ** 2
    distances = list(map(distance, zip(x_list, circle(x_list), y_list, circle(y_list))))
    i = 0
    for i in range(len(x_list)):
        if (distances[i] > threshold):
            break
    else:
        return (x_list, y_list)
    return ((x_list[0:i+1], y_list[0:i+1]), (x_list[i+1:], y_list[i+1:]), 'cut')


# Cut the x_list, y_list in TWO parts if possible to cut.
def cut_r_theta_lists(r_list, theta_list):
    threshold = 0.5
    distance = lambda l: (l[0] + l[1]) / 2 * abs(l[2] - l[3])
    distances = list(map(distance, zip(r_list, circle(r_list), theta_list, circle(theta_list))))
    i = 0
    for i in range(len(r_list)):
        if (distances[i] > threshold):
            break
    else:
        return (r_list, theta_list)
    return ((r_list[0:i+1], theta_list[0:i+1]), (r_list[i+1:], theta_list[i+1:]), 'cut')


if __name__ == '__main__':
    print('geometer.py test OK. ') if example_self_test() else print('Error met in geometer. ')
