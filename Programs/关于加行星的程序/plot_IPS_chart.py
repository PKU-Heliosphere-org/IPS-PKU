#
# plot_IPS_chart.py
# Plot IPS chart.
# Coded by Dr. L. Zhang <lifeaint1987@hotmail.com>.

from numpy import cos, sin, pi, sqrt, tan, arange
import matplotlib.pyplot as plt
from geometer import rewind_point_list, compute_cme_circle, cut_x_y_lists, cut_r_theta_lists
from astropy.coordinates import SkyCoord, GCRS, get_body, solar_system_ephemeris
from celestial_body_location import Grid_Sky
from views import *


def plot_up_head_view(list_stars_alt_az, planet_pos=None, grid=None, title=None, filename='IPS_uphead.png'):
    from matplotlib.text import Annotation
    r_plot = 1
    list_stars_seen = list(filter(lambda entry: entry[0] >= 0, list_stars_alt_az))
    list_r_theta = map(lambda entry_alt_az: cal_up_head_view_a_star(entry_alt_az[0], entry_alt_az[1], r_plot), list_stars_seen)
    list_r = [entry[0] for entry in list_r_theta]
    list_theta = [entry[1] for entry in list_r_theta]
    plt.clf()
    ax = plt.subplot(111, projection='polar')
    ax.plot(list_theta, list_r, 'b*', markersize=8)
    # ax.grid(False)
    ax.grid(True, linestyle='-')
    grid_y = cos(arange(0.0, 90.0 / 180.0 * pi, 10.0 / 180.0 * pi))
    grid_x = arange(0.0, 2 * pi, 2 * pi / 12)
    plt.xticks(grid_x, [])
    plt.yticks(grid_y, [])
    plt.ylim((0, r_plot))
    if planet_pos:
        for i_planet, planet_loc in enumerate(planet_pos): 
            sun_r, sun_theta = cal_up_head_view_a_star(planet_loc[0], planet_loc[1], r_plot)
            if (not i_planet == 0): 
                char = ['', 'Me', 'V', 'Ma', 'J', 'S'][i_planet]
                ax.text(sun_theta, sun_r, char, color='r', size=16)
            if (i_planet == 0): 
                ax.plot(sun_theta, sun_r, 'ro', markersize=20)
                for r_cme in [30, 60]:
                    loc_cme_az_alt = compute_cme_circle(az_alt_deg_sun_viewed=(planet_loc[1], planet_loc[0]), r_cme=r_cme)
                    cme_x = []
                    cme_y = []
                    for (az, alt) in loc_cme_az_alt:
                        cme_one_x, cme_one_y = cal_up_head_view_a_star(alt, az, r_plot)
                        cme_x.append(cme_one_x)
                        cme_y.append(cme_one_y)
                    ax.plot(cme_y, cme_x, 'g') # Theta must be ahead.
    if grid:
        conv_lat_lines, conv_lon_lines, text_lat_lines, text_lon_lines = grid
        for text, line in zip(text_lat_lines, conv_lat_lines):
            if not ((text == []) or (line == [])):
                wrapped_line_x, wrapped_line_y = (line[0], line[1])
                cut_line = cut_r_theta_lists(wrapped_line_x, wrapped_line_y)
                if (len(cut_line) == 2):
                    ax.plot(wrapped_line_y, wrapped_line_x, 'm--', linewidth=1)
                else:
                    (wrapped_line_x1, wrapped_line_y1), (wrapped_line_x2, wrapped_line_y2), _ = cut_line
                    ax.plot(wrapped_line_y1, wrapped_line_x1, 'm--', linewidth=1)
                    ax.plot(wrapped_line_y2, wrapped_line_x2, 'm--', linewidth=1)
        for text, line in zip(text_lon_lines, conv_lon_lines):
            if not ((text == []) or (line == [])):
                wrapped_line_x, wrapped_line_y = (line[0], line[1])
                cut_line = cut_r_theta_lists(wrapped_line_x, wrapped_line_y)
                if (len(cut_line) == 2):
                    ax.plot(wrapped_line_y, wrapped_line_x, 'k--', linewidth=1)
                else:
                    (wrapped_line_x1, wrapped_line_y1), (wrapped_line_x2, wrapped_line_y2), _ = cut_line
                    ax.plot(wrapped_line_y1, wrapped_line_x1, 'k--', linewidth=1)
                    ax.plot(wrapped_line_y2, wrapped_line_x2, 'k--', linewidth=1)
    ax.set_title(title if title else 'Sky Chart Of IPS Sources')
    plt.savefig(filename)


def plot_south_window_view(list_stars_alt_az, sun_pos=None, grid=None, title=None, filename='IPS_south_window.png'):
    from matplotlib.text import Annotation
    r_plot = 1
    list_stars_seen = list(filter(lambda entry: entry[0] >= 0, list_stars_alt_az))
    list_xy = map(lambda entry_alt_az: cal_south_window_view_a_star(entry_alt_az[0], entry_alt_az[1]), list_stars_seen)
    list_x = [entry[0] for entry in list_xy]
    list_y = [entry[1] for entry in list_xy]
    plt.clf()
    ax = plt.subplot(111)
    ax.plot(list_x, list_y, 'b*', markersize=8)
    # ax.grid(False)
    ax.grid(True, linestyle='-')
    grid_x = sin(arange(-80.0 / 180.0 * pi, 80.0 / 180.0 * pi, 10.0 / 180.0 * pi)) / sin(80.0 / 180.0 * pi)
    grid_y = tan(arange(0.0, 80.0 / 180.0 * pi, 10.0 / 180.0 * pi)) / tan(80.0 / 180.0 * pi)
    plt.xticks(grid_x, [])
    plt.yticks(grid_y, [])
    plt.xlim((-1, 1))
    plt.ylim((0, 1))
    ax.set_aspect('equal')
    if sun_pos:
        sun_x, sun_y = cal_south_window_view_a_star(sun_pos[0], sun_pos[1])
        ax.plot(sun_x, sun_y, 'ro', markersize=15)
        for r_cme in [30, 60]:
            loc_cme_az_alt = compute_cme_circle(az_alt_deg_sun_viewed=(sun_pos[1], sun_pos[0]), r_cme=r_cme)
            cme_x = []
            cme_y = []
            for (az, alt) in loc_cme_az_alt:
                cme_one_x, cme_one_y = cal_south_window_view_a_star(alt, az)
                cme_x.append(cme_one_x)
                cme_y.append(cme_one_y)
            ax.plot(cme_x, cme_y, 'g')
    if grid:
        conv_lat_lines, conv_lon_lines, text_lat_lines, text_lon_lines = grid
        for text, line in zip(text_lat_lines, conv_lat_lines):
            if not ((text == []) or (line == [])):
                ax.plot(line[0], line[1], 'm--', linewidth=1)
        for text, line in zip(text_lon_lines, conv_lon_lines):
            if not ((text == []) or (line == [])):
                ax.plot(line[0], line[1], 'k--', linewidth=1)
    ax.set_title(title if title else 'Sky Chart Of IPS Sources')
    plt.savefig(filename)


def plot_Hammer_view(list_stars_alt_az, sun_pos=None, grid=None, title=None, filename='IPS_Hammer.png'):
    from matplotlib.text import Annotation
    import operator
    list_stars_seen = list(filter(lambda entry: entry[0] >= 0, list_stars_alt_az))
    list_xy = map(lambda entry_alt_az: cal_Hammer_projection(entry_alt_az[0], entry_alt_az[1]), list_stars_seen)
    list_x = [entry[0] for entry in list_xy]
    list_y = [entry[1] for entry in list_xy]
    plt.clf()
    ax = plt.subplot(111)
    ax.plot(list_x, list_y, 'b*', markersize=8)
    ax.grid(False)
    #ax.grid(True)
    #grid_x = sin(arange(-80.0 / 180.0 * pi, 80.0 / 180.0 * pi, 10.0 / 180.0 * pi)) / sin(80.0 / 180.0 * pi)
    #grid_y = tan(arange(0.0, 80.0 / 180.0 * pi, 10.0 / 180.0 * pi)) / tan(80.0 / 180.0 * pi)
    #plt.xticks(grid_x, [])
    #plt.yticks(grid_y, [])
    plt.xlim((-2 * sqrt(2), 2 * sqrt(2)))
    plt.ylim((-1 * sqrt(2), 1 * sqrt(2)))
    ax.set_aspect('equal')
    if sun_pos:
        sun_x, sun_y = cal_Hammer_projection(sun_pos[0], sun_pos[1])
        ax.plot(sun_x, sun_y, 'ro', markersize=20)
        for r_cme in [30, 60]:
            loc_cme_az_alt = compute_cme_circle(az_alt_deg_sun_viewed=(sun_pos[1], sun_pos[0]), r_cme=r_cme)
            cme_x = []
            cme_y = []
            for (az, alt) in loc_cme_az_alt:
                cme_one_x, cme_one_y = cal_Hammer_projection(alt, az)
                cme_x.append(cme_one_x)
                cme_y.append(cme_one_y)
            ax.plot(cme_x, cme_y, 'g')
    if grid:
        conv_lat_lines, conv_lon_lines, text_lat_lines, text_lon_lines = grid
        for text, line in zip(text_lat_lines, conv_lat_lines):
            if not ((text == []) or (line == [])):
                wrapped_line_x, wrapped_line_y = rewind_point_list(line[0], line[1])
                cut_line = cut_x_y_lists(wrapped_line_x, wrapped_line_y)
                if (len(cut_line) == 2):
                    ax.plot(wrapped_line_x, wrapped_line_y, 'm--', linewidth=1)
                else:
                    (wrapped_line_x1, wrapped_line_y1), (wrapped_line_x2, wrapped_line_y2), _ = cut_line
                    ax.plot(wrapped_line_x1, wrapped_line_y1, 'm--', linewidth=1)
                    ax.plot(wrapped_line_x2, wrapped_line_y2, 'm--', linewidth=1)
        for text, line in zip(text_lon_lines, conv_lon_lines):
            if not ((text == []) or (line == [])):
                wrapped_line_x, wrapped_line_y = rewind_point_list(line[0], line[1])
                cut_line = cut_x_y_lists(wrapped_line_x, wrapped_line_y)
                if (len(cut_line) == 2):
                    ax.plot(wrapped_line_x, wrapped_line_y, 'k--', linewidth=1)
                else:
                    (wrapped_line_x1, wrapped_line_y1), (wrapped_line_x2, wrapped_line_y2), _ = cut_line
                    ax.plot(wrapped_line_x1, wrapped_line_y1, 'k--', linewidth=1)
                    ax.plot(wrapped_line_x2, wrapped_line_y2, 'k--', linewidth=1)
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
    plt.savefig(filename)


def example_self_test():
    from IPS_coord import string_to_IPS_sources
    from astropy.coordinates import EarthLocation, get_sun
    from astropy.time import Time
    import astropy.units as unit
    from celestial_body_location import current_alt_az
    font = {
            'family': 'serif', 
            'serif': 'Apple Garamond', 
            'size': 12
            }
    plt.rc('font', **font)
    the_string = '3C2 0019-00 3C12 3C26 0056-00 0106+01 0115-01 0116+082 0116+31 3C43 0128+03 3C48 3C49 0155-109 0202+15 3C67 0223+341 3C71 0258+350 CTA21 3C84 0320+05 0333+32 0347+05 0355+50 0403-132 0411+05 0428+20 3C119 3C120 3C138 0531+194 3C147 0548+16 0552+12 3C152 0605-08 3C158 0622+147 3C161 0713-024 0723+10 3C181 0732+33 3C186 0741-06 3C190 3C196 0817+183 0821+394 0834-19 3C208 3C216 3C222 3C225 0941-08 3C230 3C236 3C237 3C241 1031-11 1055+01 1055+20 1116+12 3C255 1117+14 1127-14 1136-13 3C263.1 3C267 1148-00 1213-17 3C273 3C275 1245-19 3C279 3C283 1323+32 3C286 1334-17 1345+12 1348-12 1355+01 1414+074 3C296 1414-03 3C298 1428-03 1434+03 1436-16 1452+16 1453-10 1508-05 1510-08 3C318 1518+04 1523+03 1524-13 1545-12 1548+05 1621-11 1623-228 1631-222 1638-025 1644-106 1650+004 1712-033 1730-13 1732-092 1748+031 1759+13 1819-096 1829+29 1830-210 1835+134 1858+171 1908-20 1915-121 1938-15 3C422 2120-102 2135-20 2203-18 3C446 2251+16 3C454.3 3C456 3C459 2318-16 2347-02 2354+14'
    the_IPS_sources = string_to_IPS_sources(the_string)
    Beijing = EarthLocation(lat=40.0 * unit.deg, lon=116.0 * unit.deg)
    Mingantu = EarthLocation(lat=42.3 * unit.deg, lon=115.0 * unit.deg)
    UT_offset = 8 * unit.hour
    hour = 12
 #   for day in range(1 * 365 + 1):
    for day in range(10):
        the_time = Time('2000-01-01 {:02d}:00'.format(hour)) - UT_offset + day * unit.day
        list_alt_az_IPS_sources = list(map(lambda IPS_source: current_alt_az(IPS_source, Mingantu, the_time), the_IPS_sources))
        the_sun = get_sun(the_time)
        sun_alt_az = current_alt_az(the_sun, Mingantu, the_time)
        planets = [sun_alt_az]
        for planet_name in ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn']: 
            with solar_system_ephemeris.set('builtin'):
                the_planet = get_body(planet_name, the_time, Beijing)
            the_planet_alt_az = current_alt_az(the_planet, Mingantu, the_time)
            planets.append(the_planet_alt_az)
        the_grid = Grid_Sky()
        the_grid_to_plot_up_head = convert_grid_to_r_theta_plots(the_grid, Mingantu, the_time)
        # the_grid_to_plot_south_window = convert_grid_to_r_theta_plots(the_grid, Beijing, the_time,
        #                                                              func_view=cal_south_window_view_a_star)
        # the_grid_to_plot_Hammer = convert_grid_to_r_theta_plots(the_grid, Beijing, the_time,
                                                                      #func_view=cal_Hammer_projection)
        plot_up_head_view(
            list_alt_az_IPS_sources,
            planet_pos=planets,
            grid=the_grid_to_plot_up_head,
            title='Uphead IPS Sources View In Ming''antu At {} {:02d}:00 LT'.format(the_time.iso[:10], hour),
            filename='output/IPS_uphead_{:03d}.png'.format(day)
        )
        #plot_south_window_view(
        #    list_alt_az_IPS_sources,
        #    sun_pos=sun_alt_az,
        #    grid=the_grid_to_plot_south_window,
        #    title='South Window IPS Sources View In Beijing At 2000-07-14 {:02d}:00 LT'.format(hour),
        #    filename='IPS_south_window_{:02d}.png'.format(hour)
        #)
        #plot_Hammer_view(
        #    list_alt_az_IPS_sources,
        #    sun_pos=sun_alt_az,
        #    grid=the_grid_to_plot_Hammer,
        #    title='Hammer IPS Sources View In Beijing At 2000-07-14 {:02d}:00 LT'.format(hour),
        #    filename='IPS_Hammer_{:02d}.png'.format(hour)
        #)
    return True


if __name__ == '__main__':
    print('celestial_body_location.py test OK. ') if example_self_test() else print('Error met in celestial_body_location. ')
