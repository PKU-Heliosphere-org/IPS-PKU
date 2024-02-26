#
# plot_IPS_chart.py
# Plot IPS chart.
# Coded by Dr. L. Zhang <lifeaint1987@hotmail.com>.

from numpy import cos, sin, pi, sqrt, tan, arange
import matplotlib.pyplot as plt
from geometer import rewind_point_list, compute_cme_circle, cut_x_y_lists, cut_r_theta_lists
from celestial_body_location import Grid_Sky


def cal_up_head_view_a_star(alt_deg, az_deg, r_plot=1):
    """
    :param az_deg: Azimuth in deg. North is 0; east is 90.
    :param alt_deg: Altitude in deg. Positive for what can be seen.
    :param r_plot: Radius of the plot.
    :return: a tuple (r, theta), where r is the polar radius, theta is polar angle in radians.
    """
    r = r_plot * cos(alt_deg * pi / 180)
    theta = pi / 2 + az_deg * pi / 180
    return (r, theta)


def cal_south_window_view_a_star(alt_deg, az_deg, alt_top_window_deg=80.0, az_left_deg=100.0, az_right_deg=260.0):
    """
    :param alt_deg: Altitude in deg. Positive for what can be seen.
    :param az_deg: Azimuth in deg. North is 0; east is 90.
    :param alt_top_window_deg: The altitude of top frame of the window.
    :param az_left_deg: Azimuth of left frame of the window.
    :param az_right_deg: Azimuth of right frame of the window.
    :return: a tuple (x, y), normalised to the range (-1, 1) x (0, 1): x = -1 -> left, 1 -> right; y = 0 -> bottom, 1 -> up.
    """
    az_centre = (az_left_deg + az_right_deg) / 2
    alpha = (az_deg - az_centre) * pi / 180
    x = sin(alpha) / sin((az_right_deg - az_centre) * pi / 180)
    y = tan(alt_deg * pi / 180) / (cos(alpha) * tan(alt_top_window_deg * pi / 180))
    return (x, y)


def cal_Hammer_projection(alt_deg, az_deg):
    alt_rad = alt_deg * pi / 180
    az_rad = (az_deg - 180) * pi / 180
    divider = sqrt(1 + cos(alt_rad) * cos(az_rad / 2))
    x = 2 * sqrt(2) * cos(alt_rad) * sin(az_rad / 2) / divider
    y = sqrt(2) * sin(alt_rad) / divider
    return (x, y)


def convert_grid_to_r_theta_plots(the_grid, loc, time, r_plot=1, func_view=cal_up_head_view_a_star):
    from astropy.coordinates import SkyCoord
    from celestial_body_location import current_alt_az
    import operator
    conv_lat_lines = []
    text_lat_lines = []
    conv_lon_lines = []
    text_lon_lines = []

    for lat_line in the_grid.lat_lines:
        r_list = []
        theta_list = []
        text_list = []
        for point in lat_line:
            point_as_star = SkyCoord(point[0], point[1], unit='deg')
            loc_point_seen = current_alt_az(point_as_star, loc, time)
            if (func_view == cal_up_head_view_a_star):
                if (loc_point_seen[0] > 0):
                    r, theta = func_view(loc_point_seen[0], loc_point_seen[1], r_plot)
                    r_list.append(r)
                    theta_list.append(theta)
                    text_list.append(str(point[1]))
            elif (func_view == cal_south_window_view_a_star):
                if (90 < loc_point_seen[1] < 270 and loc_point_seen[0] > 0):
                    x, y = func_view(loc_point_seen[0], loc_point_seen[1])
                    r_list.append(x)
                    theta_list.append(y)
                    text_list.append(str(point[1]))
            else:
                if (0 <= loc_point_seen[1] <= 360):
                    x, y = func_view(loc_point_seen[0], loc_point_seen[1])
                    r_list.append(x)
                    theta_list.append(y)
                    text_list.append(str(point[1]))
        if (not r_list == []):
            conv_lat_lines.append((r_list, theta_list))
            text_lat_lines.append(text_list)
    for lon_line in the_grid.lon_lines:
        r_list = []
        theta_list = []
        text_list = []
        for point in lon_line:
            point_as_star = SkyCoord(point[0], point[1], unit='deg')
            loc_point_seen = current_alt_az(point_as_star, loc, time)
            if (func_view == cal_up_head_view_a_star):
                if (loc_point_seen[0] > 0):
                    r, theta = func_view(loc_point_seen[0], loc_point_seen[1], r_plot)
                    r_list.append(r)
                    theta_list.append(theta)
                    text_list.append(str(point[0]))
            elif (func_view == cal_south_window_view_a_star):
                if (90 < loc_point_seen[1] < 270 and loc_point_seen[0] >= 0):
                    x, y = func_view(loc_point_seen[0], loc_point_seen[1])
                    r_list.append(x)
                    theta_list.append(y)
                    text_list.append(str(point[0]))
            else:
                if (0 <= loc_point_seen[1] <= 360):
                    x, y = func_view(loc_point_seen[0], loc_point_seen[1])
                    r_list.append(x)
                    theta_list.append(y)
                    text_list.append(str(point[0]))
        if (not r_list == []):
            conv_lon_lines.append((r_list, theta_list))
            text_lon_lines.append(text_list)
    return (conv_lat_lines, conv_lon_lines, text_lat_lines, text_lon_lines)


