#
# IPS_coord.py
# Convert IPS_coord string to a list of celestial bodies.
# Coded by Dr. L. Zhang <lifeaint1987@hotmail.com>.

from astropy.coordinates import SkyCoord, GCRS
from astropy.time import Time


def IPS_lon_lat_source(word):
    try:
        hour = int(word[0:2])
        minute = int(word[2:4])
        sign = 1 if word[4] == '+' else -1
        abs_lat = int(word[5:7])
        dec = int(word[7]) if len(word) > 7 else 0
        lon = (hour / 24.0 + minute / 1440.0) * 360
        lat = sign * (abs_lat + dec / 10.0)
        return SkyCoord(lon, lat, unit='deg')
    except:
        print('in lon lat')
        print(word)
        exit(-1)


def IPS_3C_source(word):
    try:
        return SkyCoord.from_name(word)
    except:
        print('in 3C')
        print(word)
        exit(-1)


def IPS_CTA_source(word):
    try:
        return SkyCoord.from_name(word)
    except:
        print('in CTA')
        print(word)
        exit(-1)


is_clean_lon_lat_word = lambda word: not (word[0:2] == '3C' or word[0:2] == 'CT')
is_clean_3C_word = lambda word: word[0:2] == '3C'
is_clean_CTA_word = lambda word: word[0:3] == 'CTA'
string_lon_lat_to_IPS_sources = lambda string: list(map(lambda w: IPS_lon_lat_source(w), filter(is_clean_lon_lat_word, string.split(' '))))
string_3C_to_IPS_sources = lambda string: list(map(lambda w: IPS_3C_source(w), filter(is_clean_3C_word, string.split(' '))))
string_CTA_to_IPS_sources = lambda string: list(map(lambda w: IPS_CTA_source(w), filter(is_clean_CTA_word, string.split(' '))))
string_to_IPS_sources = lambda string: string_lon_lat_to_IPS_sources(string) + string_3C_to_IPS_sources(string) + string_CTA_to_IPS_sources(string)


def IRS_sources(time='J2000'):
    the_string = '3C2 0019-00 3C12 3C26 0056-00 0106+01 0115-01 0116+082 0116+31 3C43 0128+03 3C48 3C49 0155-109 0202+15 3C67 0223+341 3C71 0258+350 CTA21 3C84 0320+05 0333+32 0347+05 0355+50 0403-132 0411+05 0428+20 3C119 3C120 3C138 0531+194 3C147 0548+16 0552+12 3C152 0605-08 3C158 0622+147 3C161 0713-024 0723+10 3C181 0732+33 3C186 0741-06 3C190 3C196 0817+183 0821+394 0834-19 3C208 3C216 3C222 3C225 0941-08 3C230 3C236 3C237 3C241 1031-11 1055+01 1055+20 1116+12 3C255 1117+14 1127-14 1136-13 3C263.1 3C267 1148-00 1213-17 3C273 3C275 1245-19 3C279 3C283 1323+32 3C286 1334-17 1345+12 1348-12 1355+01 1414+074 3C296 1414-03 3C298 1428-03 1434+03 1436-16 1452+16 1453-10 1508-05 1510-08 3C318 1518+04 1523+03 1524-13 1545-12 1548+05 1621-11 1623-228 1631-222 1638-025 1644-106 1650+004 1712-033 1730-13 1732-092 1748+031 1759+13 1819-096 1829+29 1830-210 1835+134 1858+171 1908-20 1915-121 1938-15 3C422 2120-102 2135-20 2203-18 3C446 2251+16 3C454.3 3C456 3C459 2318-16 2347-02 2354+14'
    the_IPS_sources = string_to_IPS_sources(the_string)
    return the_IPS_sources


def example_self_test():
    the_string = '3C2 0019-00 3C12 3C26 0056-00 0106+01 0115-01 0116+082 0116+31 3C43 0128+03 3C48 3C49 0155-109 0202+15 3C67 0223+341 3C71 0258+350 CTA21 3C84 0320+05 0333+32 0347+05 0355+50 0403-132 0411+05 0428+20 3C119 3C120 3C138 0531+194 3C147 0548+16 0552+12 3C152 0605-08 3C158 0622+147 3C161 0713-024 0723+10 3C181 0732+33 3C186 0741-06 3C190 3C196 0817+183 0821+394 0834-19 3C208 3C216 3C222 3C225 0941-08 3C230 3C236 3C237 3C241 1031-11 1055+01 1055+20 1116+12 3C255 1117+14 1127-14 1136-13 3C263.1 3C267 1148-00 1213-17 3C273 3C275 1245-19 3C279 3C283 1323+32 3C286 1334-17 1345+12 1348-12 1355+01 1414+074 3C296 1414-03 3C298 1428-03 1434+03 1436-16 1452+16 1453-10 1508-05 1510-08 3C318 1518+04 1523+03 1524-13 1545-12 1548+05 1621-11 1623-228 1631-222 1638-025 1644-106 1650+004 1712-033 1730-13 1732-092 1748+031 1759+13 1819-096 1829+29 1830-210 1835+134 1858+171 1908-20 1915-121 1938-15 3C422 2120-102 2135-20 2203-18 3C446 2251+16 3C454.3 3C456 3C459 2318-16 2347-02 2354+14'
    the_IPS_sources = string_to_IPS_sources(the_string)
    print(the_IPS_sources)
    return True


if __name__ == '__main__':
    print('celestial_body_location.py test OK. ') if example_self_test() else print('Error met in celestial_body_location. ')
