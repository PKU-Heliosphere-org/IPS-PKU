from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body, get_moon

t = Time("2014-09-22 23:22")
loc = EarthLocation.of_site('greenwich')
with solar_system_ephemeris.set('builtin'): mars = get_body('mars', t, loc)
mars_string = mars.to_string()
print(mars_string)