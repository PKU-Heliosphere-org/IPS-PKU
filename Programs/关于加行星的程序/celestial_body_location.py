#
# celestial_body_location.py
# Compute the location of a celestial body.
# Based on the example in Astropy documentation by E. Tollerud, and K. Cruz,
# at http://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html.
# Coded by Dr. L. Zhang <lifeaint1987@hotmail.com>.

from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
import astropy.units as unit


dms_to_deg = lambda dms: dms.d + dms.m / 60.0 + dms.s / 3600.0


def current_alt_az(body, loc, time):
    """
    :param body: Astropy celestial body object.
    :param loc: Astropy EarthLocation object.
    :param time: Astropy Time object.
    :return: a tuple (alt, az) in degrees.
    """
    alt_az = body.transform_to(AltAz(obstime=time, location=loc))
    return (dms_to_deg(alt_az.alt.dms), dms_to_deg(alt_az.az.dms))


class Grid_Sky:
    def generate_lon_lines(self):
        self.lon_lines = []
        for i in range(self.nx):
            res = []
            for j in range(self.ny):
                res.append((0.0 + i * self.dx, -80 + j * self.dy))
            self.lon_lines.append(res)

    def generate_lat_lines(self):
        self.lat_lines = []
        for j in range(self.ny):
            res = []
            for i in range(self.nx):
                res.append((0.0 + i * self.dx, -80 + j * self.dy))
            self.lat_lines.append(res)

    def __init__(self, dx=10, dy=10):
        self.nx = int(360.0 / dx + 0.5) + 1
        self.ny = int(160.0 / dy + 0.5) + 1
        self.dx = 360.0 / (self.nx - 1)
        self.dy = 160.0 / (self.ny - 1)
        self.generate_lon_lines()
        self.generate_lat_lines()


def example_self_test():
    from astropy.coordinates import get_sun, EarthLocation
    Beijing = EarthLocation(lat=40.0 * unit.deg, lon=116.0 * unit.deg)
    UT_offset = 8 * unit.hour
    the_time = Time('2016-12-21 12:00') - UT_offset
    the_sun = get_sun(the_time)
    sun_alt_az = current_alt_az(the_sun, Beijing, the_time)
    print('sun_alt_az: alt = {:.3f}, az = {:.3f}'.format(sun_alt_az[0], sun_alt_az[1]))
    return True


if __name__ == '__main__':
    print('celestial_body_location.py test OK. ') if example_self_test() else print('Error met in celestial_body_location. ')
