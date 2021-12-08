import numpy as np
import math
import json
from astropy import units as u
from astropy import time
from astropy.coordinates import SkyCoord
import astropy.coordinates as coordinates
from sunpy.coordinates import frames
import sunpy.coordinates
from numpy import linalg as la

# 导入数据

# Pulsar_file = open('pulsar_data_2019_0309_V2.json', 'r')
# Pulsar = eval(Pulsar_file.read())
# PSP_file = open('icrs_psp_ephemeris_data_2019_0601_1201.json', 'r')
# PSP = eval(PSP_file.read())
f1 = open("pulsar_data_2019_0309_V3.json")
# f1 = open("quasar_CTA21.json")
Pulsar = json.load(f1)
f1.close()

f2 = open("modified_psp_ICRS_ephemeris_data_2021_0101_2022_0101.json")
PSP_icrs = json.load(f2)
f2.close()

kpc2au = 3600 * 180 / math.pi * 1e3


# 列出给定时间序号time_number下，与PSP角距离小于sill的脉冲星，返回包含几个Puslar的列表
def find(time_number, still):
    # PSP_RA = PSP[time_number]['ra']
    # PSP_Dec = PSP[time_number]['dec']
    # ans储存符合条件脉冲星的数据
    ans = []
    Pulsar_num = len(Pulsar)
    # 依次读取各个脉冲星然后判断是否满足条件
    for i in range(Pulsar_num):
        '''
        if (Pulsar_RA-PSP_RA)**2 + (Pulsar_Dec-PSP_Dec)**2 <= float(sill)**2:
            ans.append(Pulsar[i])
        '''
        # elongation_PSP_Pulsar = get_elongation_PSP_Pulsar(PSP_RA, PSP_Dec, Pulsar_RA, Pulsar_Dec)
        elongation_PSP_Pulsar = get_angle_psp_earth_source(time_number, i)
        if elongation_PSP_Pulsar <= float(still):
            ans.append(Pulsar[i])
    Pulsar_close2_PSP = ans
    return Pulsar_close2_PSP


def get_elongation_PSP_Pulsar(PSP_RA, PSP_Dec, Pulsar_RA, Pulsar_Dec):
    vector_PSP = np.array([math.cos(math.radians(PSP_RA)) * math.cos(math.radians(PSP_Dec)),
                           math.sin(math.radians(PSP_RA)) * math.cos(math.radians(PSP_Dec)),
                           math.sin(math.radians(PSP_Dec))])
    vector_Pulsar = np.array([math.cos(math.radians(Pulsar_RA)) * math.cos(math.radians(Pulsar_Dec)),
                              math.sin(math.radians(Pulsar_RA)) * math.cos(math.radians(Pulsar_Dec)),
                              math.sin(math.radians(Pulsar_Dec))])
    elongation = math.degrees(math.acos(vector_PSP.dot(vector_Pulsar)))
    elongation_PSP_Pulsar = elongation
    return elongation_PSP_Pulsar


def date_english2num(date_in_english):
    month_english_num_dict = {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04", "May": "05", "Jun": "06", "Jul": "07",
                              "Aug": "08", "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"}
    return date_in_english[:5] + month_english_num_dict[date_in_english[5:8]] + date_in_english[8:]


def get_angle_psp_earth_source(time_sn, pulsar_i):
    time_str = date_english2num(PSP_icrs[time_sn]["date"]) + " " + PSP_icrs[time_sn]["time"]
    earth_coord = sunpy.coordinates.get_earth(time_str)
    earth_coord = earth_coord.transform_to("icrs")
    earth_coord.representation_type = "cartesian"
    earth_coord = np.array([earth_coord.x.value, earth_coord.y.value, earth_coord.z.value])  # 地球坐标（icrs系，纯数值，单位为AU

    psp_coord = np.array([PSP_icrs[time_sn]["x"], PSP_icrs[time_sn]["y"], PSP_icrs[time_sn]["z"]])

    source_coord = np.array([Pulsar[pulsar_i]["x_icrs"], Pulsar[pulsar_i]["y_icrs"], Pulsar[pulsar_i]["z_icrs"]])

    vec_earth_to_psp = psp_coord - earth_coord
    vec_earth_to_source = source_coord - earth_coord
    angle_psp_earth_source = math.acos(np.dot(vec_earth_to_psp, vec_earth_to_source) / la.norm(vec_earth_to_psp, 2) /
                                       la.norm(vec_earth_to_source, 2)) * 180 / math.pi
    # 地球与PSP连线 和 地球与射电源连线 之间的夹角（degree
    return angle_psp_earth_source


if __name__ == '__main__':
    # print('CTA21-Earth-PSP angle in ICRS, 20210430-0502 every 12 hours.')
    # timesn_lst = [int(17136 + i) for i in range(0, 432, 72)]
    # for timesn in timesn_lst:
    #     print(get_angle_psp_earth_source(timesn, 0))

    threshold_elongation_PSP_Pulsar = 20.0
    # the threshold of elongation between PSP and Pulsar,
    # below which the corresponding Pulsar information will be corrected
    num_times = len(PSP_icrs)  # time interval between adjacent times is 10 minutes
    print(num_times)
    PSP_and_Puslar_nearby = []
    # for i_time in range(1,num_times,6):
    # for i_time in [13513, 13393, 13537, 13249, 13105, 8809]:
    # for i_time in range(13513 - 6000, 13513 - 2400, 6):
    # for i_time in range(22752, 22896, 1):
    # for center_timesn in (22802, 38936):
    # for center_timesn in [2016]:
    for center_timesn in [17280]:
        # for i_time in range(-20, 22, 2):
        # for i_time in range(0, 720, 6):
        for i_time in range(0, 13248, 72):
            # time_number = int(center_timesn + 144 * i_time)
            time_number = int(center_timesn + i_time)
            still = threshold_elongation_PSP_Pulsar
            Pulsar_close2_PSP = find(time_number, still)
            '''
            Following lists the information of Pulsar, which is close to PSP (at a certain time number) in terms of RA & Dec coordinates
            [{'number': 405, 'name': 'B1010-23', 'jname': 'J1012-2337', 'RA': 306.2808333333333, 'Dec': -22.360444444444443, 'S400': 4.1, 'Dist': 0.98, 'Dist_date': 1978.0},
             {'number': 418, 'name': 'B1016-16', 'jname': 'J1018-1642', 'RA': 309.3358333333333, 'Dec': -15.297194444444445, 'S400': 5.1, 'S1400': 0.7, 'Dist': 25.0, 'Dist_date': 1985.0},
             {'number': 440, 'name': 'B1039-19', 'jname': 'J1041-1942', 'RA': 320.80083333333334, 'Dec': -18.296222222222223, 'pRA': -1.0, 'pDec': 14.0, 'S400': 28.0, 'S1400': 2.3, 'Dist': 2.53, 'Dist_date': 1978.0}]
            '''
            time = PSP_icrs[time_number]['date'] + ' ' + PSP_icrs[time_number]['time']
            PSP_Pulsar_now = {'time': time, 'PSP': PSP_icrs[time_number], 'Pulsar': Pulsar_close2_PSP}
            PSP_and_Puslar_nearby.append(PSP_Pulsar_now)
            print(i_time)
            '''
            Following lists the information of PSP at a certain time number
            {'time_series_number': 20000, 'jd_day': 2458774.388888889, 'date': '2019-Oct-17', 'time': '21:20:00.0000', 'ra': 311.7440386867768, 'dec': -20.12758857456606, 'distance': 0.8223222773823001}
            '''

            '''
            The function yet to be coded is to add 'Pulsar_close2_PSP' to 'PSP' as a part of the list object.
            '''
    print('down')
    jStr = json.dumps(PSP_and_Puslar_nearby, indent=2)
    # with open('PSP_and_Pulsar_nearby_20_days_around_20200607_and_0927.json', 'w') as f:
    # with open('PSP_and_Pulsar_nearby_from_2021_0115_to_0119.json', 'w') as f:
    with open('PSP_and_Pulsar_nearby_from_2021_0501_to_0731.json', 'w') as f:
        f.write(jStr)
    # file_name = 'PSP_and_Pulsar_nearby_2019_0309_modified_by_cuibo_V5.json'
    # with open(file_name, 'w') as file_object:
    #     json.dump(PSP_and_Puslar_nearby, file_object, indent=2)
