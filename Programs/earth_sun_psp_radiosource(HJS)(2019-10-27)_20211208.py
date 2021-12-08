import math
import mayavi.mlab as mlab
import json
from astropy import units as u
from astropy import time
from astropy.coordinates import SkyCoord
import astropy.coordinates as coordinates
import numpy as np
from sunpy.coordinates import frames
import sunpy.coordinates
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as plot3d
import colorsys
import julian
import datetime
from numpy import linalg as la
import xlsxwriter

# f1 = open("modified_psp_ephemeris_data_2020_0101_1231_heliocentrictrueecliptic的副本.json")
# psp_time_pos = json.load(f1)
# f1.close()

# f1 = open("modified_psp_ephemeris_data_2020_0101_1231_heliocentrictrueecliptic的副本.json")
f1 = open("modified_psp_ICRS_ephemeris_data_2021_0101_2022_0101.json")
psp_time_pos = json.load(f1)
f1.close()

# f2 = open("PSP_and_Pulsar_nearby_2019_0309(pjy)(2019-10-27).json")
# psp_radio_source_time_pos_icrs = json.load(f2)
# f2.close()

# f2 = open("PSP_and_Pulsar_nearby_30_days_around_20200927.json")
# psp_radio_source_time_pos_icrs = json.load(f2)
# f2.close()

# f2 = open("PSP_and_Pulsar_nearby_20_days_around_20200607_and_0927.json")
# psp_radio_source_time_pos_icrs = json.load(f2)
# f2.close()

# f2 = open("PSP_and_Pulsar_nearby_from_2021_0115_to_0119_from_2021_0427_to_0501.json")
# psp_radio_source_time_pos_icrs = json.load(f2)
# f2.close()

# f2 = open("PSP_and_Pulsar_3C286_nearby_from_2020_0926_to_0928.json")
f2 = open("PSP_and_Quasar_CTA21_nearby_from_2021_0430_to_0502.json")
psp_radio_source_time_pos_icrs = json.load(f2)
f2.close()

f4 = open("modified_Earth_ICRS_ephemeris_data_2021_0101_2022_0101.json")
earth_time_pos = json.load(f4)
f4.close()

f5 = open("modified_Sun_ICRS_ephemeris_data_2021_0101_2022_0101.json")
sun_time_pos = json.load(f5)
f5.close()

aa = 1*u.kpc
kpc2au = aa.to(u.AU).value
# kpc2au = 3600 * 180 / math.pi * 1e3


def date_english2num(date_in_english):
    month_english_num_dict = {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04", "May": "05", "Jun": "06", "Jul": "07",
                              "Aug": "08", "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"}
    return date_in_english[:5] + month_english_num_dict[date_in_english[5:8]] + date_in_english[8:]


# print(date_english2num("2019-Aug-31"))

"""
print(psp_time_pos[13104])
coord1 = SkyCoord(psp_time_pos[13105]["x"]*u.AU, psp_time_pos[13105]["y"]*u.AU, psp_time_pos[13105]["z"]*u.AU,
                  representation_type='cartesian', frame="heliocentrictrueecliptic", obstime="2019-08-31 00:10:00.0000")
coord2 = coord1._to(frame="icrs")
print(coord1)
print(coord2)

print(psp_radio_source_time_pos_icrs[13105//6])
print(100*"*")
coord3 = SkyCoord([116.95027976790382, 133]*u.deg, [22.455713272783722, 40]*u.deg, [0.19281336903572055, 0.35]*u.AU,
                  frame="icrs", obstime="2019-08-31 00:10:00.0000")
coord4 = coord3.transform_to("heliocentrictrueecliptic")
coord4.representation_type = "cartesian"
print(coord3)
print(coord4)
print(coord4.x.value)
print(coord4.y.value)
print(coord4.z.value)
"""

"""
earth_coord = sunpy.coordinates.get_earth("2019-08-31 00:10:00.0000")
print(earth_coord)
earth_coord = earth_coord.transform_to("heliocentrictrueecliptic")
earth_coord.representation_type = "cartesian"
print(earth_coord)
earth_coord = [earth_coord.x.value, earth_coord.y.value, earth_coord.z.value]
print(earth_coord)
print(type(earth_coord))
"""


def cal_pos_for_the_time(time_sn):
    earth_coord = np.array([earth_time_pos[time_sn]["x"], earth_time_pos[time_sn]["y"], earth_time_pos[time_sn]["z"]])  # 地球坐标（ICRS，纯数值，单位为AU

    nearby_file_index = 10000

    for i in range(len(psp_radio_source_time_pos_icrs)):
        if psp_radio_source_time_pos_icrs[i]["PSP"]["time_series_number"] == time_sn:
            nearby_file_index = i

    psp_present_coord = np.array([psp_time_pos[time_sn]["x"], psp_time_pos[time_sn]["y"], psp_time_pos[time_sn]["z"]])

    radio_source_name_lst = []
    radio_source_x_lst = []
    radio_source_y_lst = []
    radio_source_z_lst = []
    for single_source in psp_radio_source_time_pos_icrs[nearby_file_index]["Pulsar"]:
        radio_source_name_lst.append(single_source["name"])
        radio_source_x_lst.append(single_source["x_icrs"])
        radio_source_y_lst.append(single_source["y_icrs"])
        radio_source_z_lst.append(single_source["z_icrs"])
    vec_earth_to_source_normalized_lst = []
    for i in range(len(radio_source_x_lst)):
        vec_earth_to_source = np.array([radio_source_x_lst[i] - earth_coord[0],
                                        radio_source_y_lst[i] - earth_coord[1],
                                        radio_source_z_lst[i] - earth_coord[2]])
        vec_earth_to_source_normalized = np.array([vec_earth_to_source[0] / la.norm(vec_earth_to_source, 2),
                                                   vec_earth_to_source[1] / la.norm(vec_earth_to_source, 2),
                                                   vec_earth_to_source[2] / la.norm(vec_earth_to_source, 2)])
        vec_earth_to_source_normalized_lst.append(vec_earth_to_source_normalized)
    return earth_coord, psp_present_coord, vec_earth_to_source_normalized_lst, len(radio_source_x_lst), \
        radio_source_name_lst


def get_angle_source_earth_psp(time_sn):
    angle_psp_earth_source_lst = []
    nearby_file_index = 10000
    for i in range(len(psp_radio_source_time_pos_icrs)):
        if psp_radio_source_time_pos_icrs[i]["PSP"]["time_series_number"] == time_sn:
            nearby_file_index = i
    for single_source in psp_radio_source_time_pos_icrs[nearby_file_index]["Pulsar"]:
        earth_coord = np.array([earth_time_pos[time_sn]["x"], earth_time_pos[time_sn]["y"],
                                earth_time_pos[time_sn]["z"]])  # 地球坐标（ICRS，纯数值，单位为AU

        psp_present_coord = np.array(
            [psp_time_pos[time_sn]["x"], psp_time_pos[time_sn]["y"], psp_time_pos[time_sn]["z"]])
        source_coord = np.array([single_source["x_icrs"], single_source["y_icrs"], single_source["z_icrs"]])
        vec_earth_to_psp = psp_present_coord - earth_coord
        vec_earth_to_source = source_coord - earth_coord
        angle_psp_earth_source = math.acos(
            np.dot(vec_earth_to_psp, vec_earth_to_source) / la.norm(vec_earth_to_psp, 2) /
            la.norm(vec_earth_to_source, 2)) * 180 / math.pi
        angle_psp_earth_source_lst.append(angle_psp_earth_source)
    return angle_psp_earth_source_lst


def angle_x_axis_to_the_point(point_x, point_y):
    if point_y >= 0:
        the_angle = math.acos(point_x / (point_x**2 + point_y**2)**0.5)
    else:
        the_angle = 2 * math.pi - math.acos(point_x / (point_x**2 + point_y**2)**0.5)
    return the_angle


def delta_angle_x_axis_to_the_points(point1_x, point1_y, point2_x, point2_y):
    return angle_x_axis_to_the_point(point1_x, point1_y) - angle_x_axis_to_the_point(point2_x, point2_y)


def cal_p_point_r_theta(time_sn):
    coords_tuple = cal_pos_for_the_time(time_sn)
    earth_coord = coords_tuple[0]
    psp_present_coord = coords_tuple[1]
    vec_earth_to_source_normalized_lst = coords_tuple[2]
    radio_source_name_lst = coords_tuple[4]
    sun_coord = np.array([sun_time_pos[time_sn]["x"], sun_time_pos[time_sn]["y"],
                          sun_time_pos[time_sn]["z"]])  # 太阳坐标（ICRS，纯数值，单位为AU
    vec_sun_to_psp = psp_present_coord - sun_coord
    dist_sun_psp = la.norm(vec_sun_to_psp, 2)
    dist_sun_p_point_lst, p_point_elongation_angle_lst, delta_r_lst, delta_theta_lst, p_point_coord_lst = \
        [], [], [], [], []
    vec_earth_to_sun = sun_coord - earth_coord
    vec_earth_to_psp = psp_present_coord - earth_coord
    for vec_earth_to_source_normalized in vec_earth_to_source_normalized_lst:
        d = np.dot(vec_earth_to_sun, vec_earth_to_source_normalized)  # distance between Earth and P-point (AU)
        p_point_coord = earth_coord + d * vec_earth_to_source_normalized  # coordinates of P-point
        p_point_coord_lst.append(p_point_coord)
        vec_sun_to_p_point = p_point_coord - sun_coord
        dist_sun_p_point = la.norm(vec_sun_to_p_point, 2)
        dist_sun_p_point_lst.append(dist_sun_p_point)
        delta_r = dist_sun_p_point - dist_sun_psp
        # delta_r: heliocentric distance of P-point - heliocentric distance of PSP
        delta_x_axis_angle = delta_angle_x_axis_to_the_points(vec_sun_to_p_point[0], vec_sun_to_p_point[1],
                                                              vec_sun_to_psp[0], vec_sun_to_psp[1])
        if -2 * math.pi <= delta_x_axis_angle <= -math.pi or 0 <= delta_x_axis_angle <= math.pi:
            delta_theta = math.acos(np.dot(vec_sun_to_p_point, vec_sun_to_psp) / la.norm(vec_sun_to_p_point, 2) /
                                    la.norm(vec_sun_to_psp, 2)) * 180 / math.pi  # degree
        else:
            delta_theta = -math.acos(np.dot(vec_sun_to_p_point, vec_sun_to_psp) / la.norm(vec_sun_to_p_point, 2) /
                                     la.norm(vec_sun_to_psp, 2)) * 180 / math.pi  # degree
        # 这里 delta_theta 考虑了正负，正（负）为从太阳与PSP的连线经过一个小于pi的角度逆（顺）时针转到太阳与P点的连线
        delta_r_lst.append(delta_r)
        delta_theta_lst.append(delta_theta)

        vec_earth_to_p_point = p_point_coord - earth_coord
        delta_x_axis_angle = delta_angle_x_axis_to_the_points(vec_earth_to_p_point[0], vec_earth_to_p_point[1],
                                                              vec_earth_to_sun[0], vec_earth_to_sun[1])
        if -2 * math.pi <= delta_x_axis_angle <= -math.pi or 0 <= delta_x_axis_angle <= math.pi:
            p_point_elongation_angle = math.acos(np.dot(vec_earth_to_p_point, vec_earth_to_sun) /
                                                 la.norm(vec_earth_to_p_point, 2) / la.norm(vec_earth_to_sun, 2)) * \
                                       180 / math.pi  # degree
        else:
            p_point_elongation_angle = -math.acos(np.dot(vec_earth_to_p_point, vec_earth_to_sun) /
                                                  la.norm(vec_earth_to_p_point, 2) / la.norm(vec_earth_to_sun, 2)) * \
                                       180 / math.pi  # degree
        p_point_elongation_angle_lst.append(p_point_elongation_angle)
        print(delta_r, delta_theta)
    delta_x_axis_angle = delta_angle_x_axis_to_the_points(vec_earth_to_psp[0], vec_earth_to_psp[1], vec_earth_to_sun[0],
                                                          vec_earth_to_sun[1])
    if -2 * math.pi <= delta_x_axis_angle <= -math.pi or 0 <= delta_x_axis_angle <= math.pi:
        psp_elongation_angle = math.acos(np.dot(vec_earth_to_psp, vec_earth_to_sun) / la.norm(vec_earth_to_psp) /
                                         la.norm(vec_earth_to_sun)) * 180 / math.pi  # degree
    else:
        psp_elongation_angle = -math.acos(np.dot(vec_earth_to_psp, vec_earth_to_sun) / la.norm(vec_earth_to_psp) /
                                          la.norm(vec_earth_to_sun)) * 180 / math.pi  # degree
    # 这里的elongation angle与前面的delta_theta一样，也有正负
    return radio_source_name_lst, dist_sun_p_point_lst, dist_sun_psp, p_point_elongation_angle_lst, \
        psp_elongation_angle, delta_r_lst, delta_theta_lst, p_point_coord_lst


def write_r_theta_etc_to_excel():
    # timesn1_lst = [int(22802 + 144 * i) for i in range(-20, 22, 2)]
    # timesn2_lst = [int(38936 + 144 * i) for i in range(-20, 22, 2)]
    timesn1_lst = [int(38736 + i) for i in range(0, 576, 6)]
    timesn2_lst = [int(2016 + i) for i in range(0, 720, 6)]
    timesn3_lst = [int(16704 + i) for i in range(0, 720, 6)]
    timesn4_lst = [int(17280 + i) for i in range(0, 13248, 72)]
    timesn5_lst = [int(17136 + i) for i in range(0, 432, 72)]
    # workbook = xlsxwriter.Workbook("PSP and RadioSource information_from_2020_0926_to_0929_from_2021_0115_to_0119_from_2021_0427_to_0501.xlsx")
    # workbook = xlsxwriter.Workbook("PSP and RadioSource information_from_2020_0926_to_0929.xlsx")
    # workbook = xlsxwriter.Workbook("PSP and RadioSource information_from_2021_0115_to_0119_from_2021_0427_to_0501.xlsx")
    # workbook = xlsxwriter.Workbook("PSP and RadioSource information_from_2021_0501_to_0731.xlsx")
    workbook = xlsxwriter.Workbook("PSP and RadioSource information_from_2021_0430_to_0502.xlsx")
    worksheet = workbook.add_worksheet()
    format1 = workbook.add_format({"font_name": "Times New Roman", "font_size": 18, "align": "center",
                                   "valign": "vcenter"})
    format2 = workbook.add_format({"font_name": "Times New Roman", "font_size": 16, "align": "center",
                                   "valign": "vcenter"})
    worksheet.write_row(0, 1, ("date and time", "Heliocentric distance of PSP\n[AU]",
                               "elongation angle of PSP\n[degree]", "Radio Source name",
                               "angle source-Earth-PSP\n[degree]", "Heliocentric distance of P-point\n[AU]",
                               "elongation angle of P-point\n[degree]", "delta r\n[AU]", "phi\n[degree]"),
                        cell_format=format1)
    worksheet.set_row(0, 60)
    worksheet.set_column(0, 9, 35)
    worksheet.set_column(2, 2, 40)
    worksheet.set_column(4, 4, 35)
    worksheet.set_column(5, 5, 45)
    worksheet.set_column(6, 6, 45)
    worksheet.set_column(8, 8, 15)
    worksheet.set_column(9, 9, 25)
    # center_tuple = ("2020-06-07\nperihelion", "2020-09-27\nperihelion")
    # center_tuple = ["from 2020-09-26 to 09-29"]
    # center_tuple = ["from 2021-01-15 to 01-19", "from 2021-04-27 to 05-01"]
    # center_tuple = ["from 2021-05-01 to 07-31"]
    center_tuple = ["from 2021-04-30 to 05-02"]
    i_center = 0
    i_row = 1
    center_first_row = i_row
    single_timesn_first_row = i_row
    # for timesn_lst in (timesn1_lst, timesn2_lst):
    # for single_timesn_lst in [timesn1_lst]:
    # for single_timesn_lst in [timesn2_lst, timesn3_lst]:
    for single_timesn_lst in [timesn5_lst]:
        for timesn in single_timesn_lst:
            print(timesn)
            coords_tuple = cal_p_point_r_theta(timesn)
            radio_source_name_lst, dist_sun_p_point_lst, dist_sun_psp, p_point_elongation_angle_lst, \
                psp_elongation_angle, delta_r_lst, delta_theta_lst = coords_tuple[:7]
            angle_source_earth_psp_lst = get_angle_source_earth_psp(timesn)
            for i in range(len(radio_source_name_lst)):
                worksheet.set_row(i_row, 30)
                worksheet.write_row(i_row, 4, (radio_source_name_lst[i], float("%.4f" % angle_source_earth_psp_lst[i]),
                                               float("%.4f" % dist_sun_p_point_lst[i]),
                                               float("%.4f" % p_point_elongation_angle_lst[i]),
                                               float("%.4f" % delta_r_lst[i]), float("%.4f" % delta_theta_lst[i])),
                                    cell_format=format2)
                i_row += 1
            single_timesn_last_row = i_row - 1
            time_str = date_english2num(psp_time_pos[timesn]["date"]) + " " + psp_time_pos[timesn]["time"]
            if single_timesn_last_row > single_timesn_first_row:
                worksheet.merge_range(single_timesn_first_row, 1, single_timesn_last_row, 1, time_str,
                                      cell_format=format2)
                worksheet.merge_range(single_timesn_first_row, 2, single_timesn_last_row, 2,
                                      float("%.4f" % dist_sun_psp),
                                      cell_format=format2)
                worksheet.merge_range(single_timesn_first_row, 3, single_timesn_last_row, 3,
                                      float("%.4f" % psp_elongation_angle),
                                      cell_format=format2)
            else:
                worksheet.write_row(single_timesn_first_row, 1, (time_str, float("%.4f" % dist_sun_psp),
                                                                 float("%.4f" % psp_elongation_angle)),
                                    cell_format=format2)
            single_timesn_first_row = i_row
        center_last_row = i_row - 1
        worksheet.merge_range(center_first_row, 0, center_last_row, 0, center_tuple[i_center],
                              cell_format=format2)
        i_center += 1
        center_first_row = i_row
    workbook.close()


def plot_earth_sun_psp_radio_source(time_series_number_lst, num_of_pos_psp):
    datetime_str_lst = []
    for time_sn in time_series_number_lst:  # time_sn 应从 nearby 文件中选
        print(len(time_series_number_lst))
        time_str = date_english2num(psp_time_pos[time_sn]["date"]) + " " + psp_time_pos[time_sn]["time"]
        print(time_str)
        datetime_str = time_str.replace(':', '')
        datetime_str = datetime_str.replace('-', '')
        datetime_str = datetime_str[:-5]
        print('datetime_str: ', datetime_str)
        datetime_str_lst.append(datetime_str)

        coords_tuple = cal_pos_for_the_time(time_sn)
        coords_tuple2 = cal_p_point_r_theta(time_sn)

        earth_coord = coords_tuple[0]

        psp_time_sn_beginning = time_sn - num_of_pos_psp
        psp_time_sn_end = time_sn + num_of_pos_psp  # PSP轨迹的起点和终点对应的time_series_number
        psp_traj_x_lst = []
        psp_traj_y_lst = []
        psp_traj_z_lst = []
        for time_sn_psp_traj in range(psp_time_sn_beginning, psp_time_sn_end, 144):  # time_sn_psp_traj 是PSP轨迹上各点
            # 对应的time_series_number
            psp_traj_x_lst.append(psp_time_pos[time_sn_psp_traj]["x"])
            psp_traj_y_lst.append(psp_time_pos[time_sn_psp_traj]["y"])
            psp_traj_z_lst.append(psp_time_pos[time_sn_psp_traj]["z"])

        dtime_psp_time_pos = 10.0  # unit: min
        dtime_psp_radio_source = 60.0  # unit: min
        # nearby_file_index = int(time_sn // (dtime_psp_radio_source/dtime_psp_time_pos))
        # time_sn 对应的 nearby 文件对应的时间的index

        psp_present_coord = coords_tuple[1]

        figure = mlab.figure('Earth_Sun_PSP_RadioSource')
        figure.scene.disable_render = True  # Super duper trick

        mlab.figure(bgcolor=(0., 0., 0.), size=(1000, 800))

        sun_coord = np.array([sun_time_pos[time_sn]["x"], sun_time_pos[time_sn]["y"],
                              sun_time_pos[time_sn]["z"]])  # 太阳坐标（ICRS，纯数值，单位为AU

        mlab.points3d(sun_coord[0], sun_coord[1], sun_coord[2], scale_factor=0.1, color=(1, 1, 0))  # 画太阳
        mlab.points3d(earth_coord[0], earth_coord[1], earth_coord[2], scale_factor=0.1, color=(0, 0, 1))  # 画地球
        mlab.points3d(psp_present_coord[0], psp_present_coord[1], psp_present_coord[2], scale_factor=0.1,
                      color=(1, 0, 0))
        for i in range(len(psp_traj_x_lst)):
            mlab.points3d(psp_traj_x_lst[i], psp_traj_y_lst[i], psp_traj_z_lst[i], scale_factor=0.05, color=(1, 1, 1))
        # mlab.points3d(psp_traj_x_lst, psp_traj_y_lst, psp_traj_z_lst, color=(1, 1, 1))  # 画PSP的轨迹（离散的点

        num_radiosources = coords_tuple[3]
        dh_of_hsv = 1. / num_radiosources
        radio_source_name_lst = coords_tuple[4]
        for i in range(num_radiosources):
            h_of_hsv_tmp = dh_of_hsv * i
            s_of_hsv_tmp = 0.4
            v_of_hsv_tmp = 1.
            color_hsv = [h_of_hsv_tmp, s_of_hsv_tmp, v_of_hsv_tmp]
            color_rgb = colorsys.hsv_to_rgb(color_hsv[0], color_hsv[1], color_hsv[2])

            """崔博修改，修改了画地球到射电源的连线的部分 20191027"""
            length_los = 1.7
            vec_earth_to_source_normalized = coords_tuple[2][i]
            p_point_coord = coords_tuple2[7][i]
            mlab.plot3d([earth_coord[0], earth_coord[0] + length_los * vec_earth_to_source_normalized[0]],
                        [earth_coord[1], earth_coord[1] + length_los * vec_earth_to_source_normalized[1]],
                        [earth_coord[2], earth_coord[2] + length_los * vec_earth_to_source_normalized[2]],
                        tube_radius=0.01, color=color_rgb)
            mlab.text3d(earth_coord[0] + length_los * vec_earth_to_source_normalized[0],
                        earth_coord[1] + length_los * vec_earth_to_source_normalized[1],
                        earth_coord[2] + length_los * vec_earth_to_source_normalized[2], radio_source_name_lst[i],
                        line_width=5, scale=(0.12 / 2, 0.12 / 2, 0.12 / 2), color=color_rgb)
            mlab.points3d(p_point_coord[0], p_point_coord[1], p_point_coord[2], scale_factor=0.049,
                          color=(0, 205/255, 102/255))
            # mlab.points3d(p_point_coord[0], p_point_coord[1], p_point_coord[2], mode="cube", scale_factor=0.055,
            #               color=(1, 1, 1))
            # 画 P-point
            """结束"""

            """
            mlab.plot3d([earth_coord[0], radio_source_x_lst[i]], [earth_coord[1], radio_source_y_lst[i]],
                        [earth_coord[2], radio_source_z_lst[i]], color=color_rgb)
            mlab.text3d(radio_source_x_lst[i], radio_source_y_lst[i], radio_source_z_lst[i], radio_source_name_lst[i],
                        scale=(0.1/2, 0.1/2, 0.1/2),color=color_rgb)
            """
        figure.scene.disable_render = False  # Super duper trick
        # dir_fig = '/E-Work/Data_Analysis/FAST_IPS_PSP_Project/Figures/'
        # dir_fig = '/Users/cuibo/Desktop/My_Files/太阳风小小组/PSP_IPS/FAST_PSP_radio_source_IPS_project_20200512/' \
        #           'Earth_Sun_PSP_RadioSources_from_2020_0604_to_0609/azi_0_ele_45/'
        dir_fig = ''
        azimuth_view = 45.0
        # elevation_view = 45.0
        elevation_view = 0.0
        # azimuth_view = 225.0
        # elevation_view = 60.0
        azimuth_view_str = str(azimuth_view)
        elevation_view_str = str(elevation_view)
        focal_point_lst = [0, 0, 0]
        file_fig = 'Earth_Sun_PSP_CTA21 ' + \
                   '(datetime=' + datetime_str + ')' + \
                   '(azi&ele=' + azimuth_view_str + ',' + elevation_view_str + ')' + \
                   '.png'
        title_str1 = 'Earth & Sun & PSP & Radio Sources'
        title_str2 = '(datetime=' + datetime_str + ')' + \
                     '(azi&ele=' + azimuth_view_str + ',' + elevation_view_str + ')'
        # mlab.title(title_str,color=(0.8,0.8,0.8),size=0.2)
        mlab.text(0.065, 0.85, title_str1, color=(0.9, 0.9, 0.9), width=0.9)
        mlab.text(0.065, 0.75, title_str2, color=(0.9, 0.9, 0.9), width=0.9)
        mlab.view(azimuth=azimuth_view, elevation=elevation_view, distance=3.9, focalpoint=focal_point_lst)
        mlab.orientation_axes()
        mlab.savefig(dir_fig + file_fig)
        # mlab.show()

        # _ = input("Press [enter] to continue.")  # wait for input from the user
        mlab.close()  # close the figure to show the next one.


if __name__ == "__main__":
    timesn_lst = [int(17136 + i) for i in range(0, 432, 72)]
    plot_earth_sun_psp_radio_source(timesn_lst, 4320)
    # write_r_theta_etc_to_excel()
    # for i in timesn_lst:
    #     print(get_angle_source_earth_psp(i))

    # print(cal_p_point_r_theta(int(22802 + 144 * 0)))
    # timesn_lst = [int(38736 + i) for i in range(0, 576, 6)] + [int(2016 + i) for i in range(0, 720, 6)] + \
    #              [int(16704 + i) for i in range(0, 720, 6)]
    # timesn_lst = [int(22320 + i) for i in range(0, 12, 12)]
    # write_r_theta_etc_to_excel()
#
#     PSP_ephemeris_filename = 'psp_ephemeris_data_2019_0601_1201_heliocentrictrueecliptic的副本.json'
#     f = open(PSP_ephemeris_filename,'r')
#     PSP_time_pos_entries = json.load(f)
#     f.close()
#     lenth = len(PSP_time_pos_entries)
#     num_eff_entries = 0
#     # for i_entry in range(0,lenth,1):
#     for i_entry in range(7513, 15308):
#         num_eff_entries = num_eff_entries + 1
#         PSP_entry = PSP_time_pos_entries[i_entry]
#         date_tmp = PSP_entry['date']
#         time_tmp = PSP_entry['time']
#         datetime_tmp = date_tmp + ' '+time_tmp
#         datetime_tmp = datetime_tmp[:-5]
#         print(datetime_tmp)
#         if num_eff_entries == 1:
#             datetime_vect = [datetime_tmp]
#             i_entry_vect = i_entry
#         if num_eff_entries > 1:
#             datetime_vect.append(datetime_tmp)
#             i_entry_vect = np.append(i_entry_vect,i_entry)
#     print(datetime_vect)
#     print('i_entry_vect: ')
#     print(i_entry_vect)
#     num_times = 0
#     for datetime_tmp in datetime_vect:
#         print('datetime_tmp',datetime_tmp)
#         datetime_tmp = datetime.datetime.strptime(datetime_tmp, '%Y-%b-%d %H:%M:%S')
#         jd_tmp = julian.to_jd(datetime_tmp, fmt='jd')
#         print(jd_tmp)
#         num_times = num_times + 1
#         if num_times == 1:
#             jd_vect = jd_tmp
#         if num_times > 1:
#             jd_vect = np.append(jd_vect,jd_tmp)
#     print(jd_vect)
#
#     # datetime_str_beg = '2019-Sep-02 04:00:00'
#     # datetime_str_end = '2019-Sep-02 05:00:00'
#     datetime_str_beg = '2019-Sep-15 07:00:00'
#     datetime_str_end = '2019-Sep-15 08:00:00'
#     datetime_beg = datetime.datetime.strptime(datetime_str_beg, '%Y-%b-%d %H:%M:%S')
#     datetime_end = datetime.datetime.strptime(datetime_str_end, '%Y-%b-%d %H:%M:%S')
#     jd_beg = julian.to_jd(datetime_beg, fmt='jd')
#     jd_end = julian.to_jd(datetime_end, fmt='jd')
#     print(jd_beg, jd_end)
#     index_vect = np.arange(len(jd_vect))
#     print(index_vect)
#     #    index_in_range = index_vect(jd_vect > jd_beg) # & jd_vect < jd_end)
#     index_in_range = np.where((jd_vect > jd_beg) & (jd_vect < jd_end))
#     print(index_in_range)
#     print(jd_vect[index_in_range])
#     i_entry_select_vect = i_entry_vect[index_in_range]
#     print(i_entry_select_vect)
#     print(i_entry_select_vect[0:1])
#     _ = input("Press [enter] to continue.")  # wait for input from the user
#     plot_earth_sun_psp_radio_source(i_entry_select_vect[0:1], 4320)
#     # plot_earth_sun_psp_radio_source([13105], 4320)  # 201908310010
#     # plot_earth_sun_psp_radio_source([13249], 4320)  # 201909010010
#     # plot_earth_sun_psp_radio_source([13393], 4320)  # 201909020010
