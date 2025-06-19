import os
from pathlib import Path


PROJECT_ROOT = Path(os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../')))
SECRETS_PATH = PROJECT_ROOT/'secrets.json'

# SRTM tiles that cover Guinea or nearby countries (for robustness)
SRTM_GIN_TILE_IDS = [
    # N12 W016 to W008
    'N12W016',
    'N12W015',
    'N12W014',
    'N12W013',
    'N12W012',
    'N12W011',
    'N12W010',
    'N12W009',
    'N12W008',

    # N11 W016 to W008
    'N11W016',
    'N11W015',
    'N11W014',
    'N11W013',
    'N11W012',
    'N11W011',
    'N11W010',
    'N11W009',
    'N11W008',

    # N10 W016 to W008
    'N10W016',
    'N10W015',
    'N10W014',
    'N10W013',
    'N10W012',
    'N10W011',
    'N10W010',
    'N10W009',
    'N10W008',

    # N09 W015 to W008
    'N09W015',
    'N09W014',
    'N09W013',
    'N09W012',
    'N09W011',
    'N09W010',
    'N09W009',
    'N09W008',

    # N08 W014 to W008
    'N08W014',
    'N08W013',
    'N08W012',
    'N08W011',
    'N08W010',
    'N08W009',
    'N08W008',
    'N04W008',
    'N04W007',

    # N07 W014 to W008
    'N07W014',
    'N07W013',
    'N07W012',
    'N07W011',
    'N07W010',
    'N07W009',

    # N06 W012 to W008
    'N06W012',
    'N06W011',
    'N06W010',
    'N06W009',
    'N06W008',

    # N05 W011 to W008
    'N05W011',
    'N05W010',
    'N05W009',
    'N05W008',

    # N04 W010 to W008
    'N04W010',
    'N04W009',
    'N04W008',
]

#: Transmitter CSV files must have these header columns
REQUIRED_TRANSMITTER_FIELDS = [
    'network_name',
    'site_name',
    'latitude',  # WGS84 float
    'longitude',  # WGS84 float
    'antenna_height',  # meters
    'polarization',  # 0 (horizontal) or 1 (vertical)
    'frequency',  # megaherz
    'power_eirp',  # watts
]

#: SPLAT! Earth dielectric constant.
#: According to the SPLAT! documentation, typical Earth dielectric constants and conductivities are:
#: Salt water, 80, 5.000;
#: Good ground, 25, 0.020;
#: Fresh water, 80, 0.010;
#: Marshy land, 12, 0.007;
#: Farmland or forest, 15, 0.005;
#: Average ground, 15, 0.005;
#: Mountain or sand, 13, 0.002;
#: City, 5, 0.001;
#: Poor ground, 4, 0.001;
EARTH_DIELECTRIC_CONSTANT = 15

#: SPLAT! Earth earth_conductivity in Siemens per meter
EARTH_CONDUCTIVITY = 0.005

#: SPLAT! radio climate codes.
#: 1=Equatorial (Congo);
#: 2=Continental Subtropical (Sudan);
#: 3=Maritime Subtropical (West coast of Africa);
#: 4=Desert (Sahara);
#: 5=Continental Temperate;
#: 6=Maritime Temperate, over land (UK and west coasts of US & EU);
#: 7=Maritime Temperate, over sea
RADIO_CLIMATE = 3

#: SPLAT! time variability parameter
FRACTION_OF_TIME = 0.5

#: SPLAT! location variability parameter
FRACTION_OF_SITUATIONS = 0.5

#: SPLAT receiver sensitivity parameter in decibel-milliwatts (dBm).
#: For example, minimum received signal power of wireless networks (802.11 variants) is -100 dBm.
RECEIVER_SENSITIVITY = -110
#: WGS84 semimajor axis in meters
WGS84_A = 6378137
#: WGS84 flattening
WGS84_F = 1/298.257223563
#: WGS84 eccentricity squared (e^2)
WGS84_E2 = 2*WGS84_F - WGS84_F**2
#: Distance in meters of a geostationary satellite from the center of the Earth (and hence the center of the WGS84 ellipsoid);
#: taken from the Wikipedia article `Geostationary orbit <https://en.wikipedia.org/wiki/Geostationary_orbit>`_
R_S = 42164000
#: Distance in meters of a geostationary satellite from the WGS84 ellipsoid
H_S = R_S - WGS84_A
