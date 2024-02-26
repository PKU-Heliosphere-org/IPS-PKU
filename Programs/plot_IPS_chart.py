#
# plot_IPS_chart.py
# Plot IPS chart.
# Coded by Dr. L. Zhang <lifeaint1987@hotmail.com>.

import numpy as np
import matplotlib.pyplot as plt
from geometer import rewind_point_list
from celestial_body_location import Grid_Sky


def cal_up_head_view_a_star(alt_deg, az_deg, r_plot=1):
    """
    :param az_deg: Azimuth in deg. North is 0; east is 90.
    :param alt_deg: Altitude in deg. Positive for what can be seen.
    :param r_plot: Radius of the plot.
    :return: a tuple (r, theta), where r is the polar radius, theta is polar angle in radians.
    """
    r = r_plot * np.cos(alt_deg * np.pi / 180)
    theta = np.pi / 2 + az_deg * np.pi / 180
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
    alpha = (az_deg - az_centre) * np.pi / 180
    x = np.sin(alpha) / np.sin((az_right_deg - az_centre) * np.pi / 180)
    y = np.tan(alt_deg * np.pi / 180) / (np.cos(alpha) * np.tan(alt_top_window_deg * np.pi / 180))
    return (x, y)


def cal_Hammer_projection(alt_deg, az_deg):
    alt_rad = alt_deg * np.pi / 180
    az_rad = (az_deg - 180) * np.pi / 180
    divider = np.sqrt(1 + np.cos(alt_rad) * np.cos(az_rad / 2))
    x = 2 * np.sqrt(2) * np.cos(alt_rad) * np.sin(az_rad / 2) / divider
    y = np.sqrt(2) * np.sin(alt_rad) / divider
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


def plot_up_head_view(list_stars_alt_az, sun_pos=None, grid=None, title=None):
    from matplotlib.text import Annotation
    r_plot = 1
    list_stars_seen = list(filter(lambda entry: entry[0] >= 0, list_stars_alt_az))
    list_r_theta = list(map(lambda entry_alt_az: cal_up_head_view_a_star(entry_alt_az[0], entry_alt_az[1], r_plot), list_stars_seen))
    list_r = [entry[0] for entry in list_r_theta]
    list_theta = [entry[1] for entry in list_r_theta]
    plt.figure()
    ax = plt.subplot(111, projection='polar')
    ax.plot(list_theta, list_r, 'b*', markersize=8)
    # ax.grid(False)
    ax.grid(True, linestyle='-')
    grid_y = np.cos(np.arange(0.0, 90.0 / 180.0 * np.pi, 10.0 / 180.0 * np.pi))
    grid_x = np.arange(0.0, 2 * np.pi, 2 * np.pi / 12)
    plt.xticks(grid_x, [])
    plt.yticks(grid_y, [])
    if sun_pos:
        sun_r, sun_theta = cal_up_head_view_a_star(sun_pos[0], sun_pos[1], r_plot)
        ax.plot(sun_theta, sun_r, 'ro', markersize=20)
    if grid:
        conv_lat_lines, conv_lon_lines, text_lat_lines, text_lon_lines = grid
        for text, line in zip(text_lat_lines, conv_lat_lines):
            if not ((text == []) or (line == [])):
                ax.plot(line[1], line[0], 'm--', linewidth=1)
        for text, line in zip(text_lon_lines, conv_lon_lines):
            if not ((text == []) or (line == [])):
                ax.plot(line[1], line[0], 'k--', linewidth=1)
    ax.set_title(title if title else 'Sky Chart Of IPS Sources')
    plt.savefig('IPS_uphead.png')


def plot_south_window_view(list_stars_alt_az, sun_pos=None, grid=None, title=None):
    from matplotlib.text import Annotation
    r_plot = 1
    list_stars_seen = list(filter(lambda entry: entry[0] >= 0, list_stars_alt_az))
    list_xy = list(map(lambda entry_alt_az: cal_south_window_view_a_star(entry_alt_az[0], entry_alt_az[1]), list_stars_seen))
    list_x = [entry[0] for entry in list_xy]
    list_y = [entry[1] for entry in list_xy]
    plt.figure()
    ax = plt.subplot(111)
    ax.plot(list_x, list_y, 'b*', markersize=8)
    # ax.grid(False)
    ax.grid(True, linestyle='-')
    grid_x = np.sin(np.arange(-80.0 / 180.0 * np.pi, 80.0 / 180.0 * np.pi, 10.0 / 180.0 * np.pi)) / np.sin(80.0 / 180.0 * np.pi)
    grid_y = np.tan(np.arange(0.0, 80.0 / 180.0 * np.pi, 10.0 / 180.0 * np.pi)) / np.tan(80.0 / 180.0 * np.pi)
    plt.xticks(grid_x, [])
    plt.yticks(grid_y, [])
    plt.xlim((-1, 1))
    plt.ylim((0, 1))
    ax.set_aspect('equal')
    if sun_pos:
        sun_x, sun_y = cal_south_window_view_a_star(sun_pos[0], sun_pos[1])
        ax.plot(sun_x, sun_y, 'ro', markersize=20)
    if grid:
        conv_lat_lines, conv_lon_lines, text_lat_lines, text_lon_lines = grid
        for text, line in zip(text_lat_lines, conv_lat_lines):
            if not ((text == []) or (line == [])):
                ax.plot(line[0], line[1], 'm--', linewidth=1)
        for text, line in zip(text_lon_lines, conv_lon_lines):
            if not ((text == []) or (line == [])):
                ax.plot(line[0], line[1], 'k--', linewidth=1)
    ax.set_title(title if title else 'Sky Chart Of IPS Sources')
    plt.savefig('IPS_south_window.png')


def plot_Hammer_view(list_stars_alt_az, sun_pos=None, grid=None, title=None):
    from matplotlib.text import Annotation
    import operator
    list_stars_seen = list(filter(lambda entry: entry[0] >= 0, list_stars_alt_az))
    list_xy = list(map(lambda entry_alt_az: cal_Hammer_projection(entry_alt_az[0], entry_alt_az[1]), list_stars_seen))
    list_x = [entry[0] for entry in list_xy]
    list_y = [entry[1] for entry in list_xy]
    plt.figure()
    ax = plt.subplot(111)
    ax.plot(list_x, list_y, 'b*', markersize=8)
    ax.grid(False)
    #ax.grid(True)
    #grid_x = np.sin(np.arange(-80.0 / 180.0 * np.pi, 80.0 / 180.0 * np.pi, 10.0 / 180.0 * np.pi)) / np.sin(80.0 / 180.0 * np.pi)
    #grid_y = np.tan(np.arange(0.0, 80.0 / 180.0 * np.pi, 10.0 / 180.0 * np.pi)) / np.tan(80.0 / 180.0 * np.pi)
    #plt.xticks(grid_x, [])
    #plt.yticks(grid_y, [])
    plt.xlim((-2 * np.sqrt(2), 2 * np.sqrt(2)))
    plt.ylim((-1 * np.sqrt(2), 1 * np.sqrt(2)))
    ax.set_aspect('equal')
    if sun_pos:
        sun_x, sun_y = cal_Hammer_projection(sun_pos[0], sun_pos[1])
        ax.plot(sun_x, sun_y, 'ro', markersize=20)
    if grid:
        conv_lat_lines, conv_lon_lines, text_lat_lines, text_lon_lines = grid
        for text, line in zip(text_lat_lines, conv_lat_lines):
            if not ((text == []) or (line == [])):
                wrapped_line_x, wrapped_line_y = rewind_point_list(line[0], line[1])
                ax.plot(wrapped_line_x, wrapped_line_y, 'm--', linewidth=1)
        for text, line in zip(text_lon_lines, conv_lon_lines):
            if not ((text == []) or (line == [])):
                wrapped_line_x, wrapped_line_y = rewind_point_list(line[0], line[1])
                ax.plot(wrapped_line_x, wrapped_line_y, 'k--', linewidth=1)
    # Plot alt / az.
    for alt in [-60, -30, 0, 30, 60]:
        az_list = range(0, 370, 10)
        line_x = []
        line_y = []
        for az in az_list:
            x, y = cal_Hammer_projection(alt, az)
            line_x.append(x)
            line_y.append(y)
        ax.plot(line_x, line_y, 'k', linewidth=1)
    for az in [0, 60, 120, 180, 240, 300, 360]:
        alt_list = range(-90, 100, 10)
        line_x = []
        line_y = []
        for alt in alt_list:
            x, y = cal_Hammer_projection(alt, az)
            line_x.append(x)
            line_y.append(y)
        ax.plot(line_x, line_y, 'k', linewidth=(3 if az == 0 or az == 360 else 1))

    ax.set_title(title if title else 'Sky Chart Of IPS Sources')
    plt.savefig('IPS_Hammer.png')


def example_self_test():
    from IPS_coord import string_to_IPS_sources
    from astropy.coordinates import EarthLocation, get_sun
    from astropy.time import Time
    import astropy.units as unit
    from celestial_body_location import current_alt_az
    the_string = '3C2 0019-00 3C12 3C26 0056-00 0106+01 0115-01 0116+082 0116+31 3C43 0128+03 3C48 3C49 0155-109 0202+15 3C67 0223+341 3C71 0258+350 CTA21 3C84 0320+05 0333+32 0347+05 0355+50 0403-132 0411+05 0428+20 3C119 3C120 3C138 0531+194 3C147 0548+16 0552+12 3C152 0605-08 3C158 0622+147 3C161 0713-024 0723+10 3C181 0732+33 3C186 0741-06 3C190 3C196 0817+183 0821+394 0834-19 3C208 3C216 3C222 3C225 0941-08 3C230 3C236 3C237 3C241 1031-11 1055+01 1055+20 1116+12 3C255 1117+14 1127-14 1136-13 3C263.1 3C267 1148-00 1213-17 3C273 3C275 1245-19 3C279 3C283 1323+32 3C286 1334-17 1345+12 1348-12 1355+01 1414+074 3C296 1414-03 3C298 1428-03 1434+03 1436-16 1452+16 1453-10 1508-05 1510-08 3C318 1518+04 1523+03 1524-13 1545-12 1548+05 1621-11 1623-228 1631-222 1638-025 1644-106 1650+004 1712-033 1730-13 1732-092 1748+031 1759+13 1819-096 1829+29 1830-210 1835+134 1858+171 1908-20 1915-121 1938-15 3C422 2120-102 2135-20 2203-18 3C446 2251+16 3C454.3 3C456 3C459 2318-16 2347-02 2354+14'
    the_IPS_sources = string_to_IPS_sources(the_string)
    Beijing = EarthLocation(lat=40.0 * unit.deg, lon=116.0 * unit.deg)
    UT_offset = 8 * unit.hour
    the_time = Time('2000-07-14 14:00') - UT_offset
    list_alt_az_IPS_sources = list(map(lambda IPS_source: current_alt_az(IPS_source, Beijing, the_time), the_IPS_sources))
    the_sun = get_sun(the_time)
    sun_alt_az = current_alt_az(the_sun, Beijing, the_time)
    the_grid = Grid_Sky()
    the_grid_to_plot_up_head = convert_grid_to_r_theta_plots(the_grid, Beijing, the_time)
    the_grid_to_plot_south_window = convert_grid_to_r_theta_plots(the_grid, Beijing, the_time,
                                                                  func_view=cal_south_window_view_a_star)
    the_grid_to_plot_Hammer = convert_grid_to_r_theta_plots(the_grid, Beijing, the_time,
                                                                  func_view=cal_Hammer_projection)
    plot_up_head_view(
        list_alt_az_IPS_sources,
        sun_pos=sun_alt_az,
        grid=the_grid_to_plot_up_head,
        title='Uphead IPS Sources View In Beijing At 2000-07-14 14:00 LT'
    )
    plot_south_window_view(
        list_alt_az_IPS_sources,
        sun_pos=sun_alt_az,
        grid=the_grid_to_plot_south_window,
        title='South Window IPS Sources View In Beijing At 2000-07-14 14:00 LT'
    )
    plot_Hammer_view(
        list_alt_az_IPS_sources,
        sun_pos=sun_alt_az,
        grid=the_grid_to_plot_Hammer,
        title='Hammer IPS Sources View In Beijing At 2000-07-14 14:00 LT'
    )
    return True


if __name__ == '__main__':
    example_self_test()
    # print('celestial_body_location.py test OK. ') if example_self_test() else print('Error met in celestial_body_location. ')
