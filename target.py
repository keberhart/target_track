#!/usr/bin/env python
#
# create a list of points so we can track a target with the new Orbit ACU
#
#   9APR19 - Kyle Eberhart
#   5MAY20 - Kyle Eberhart - extend this to other targets and antennas
#
#   This project made use of Skyfield, http://rhodesmill.org/skyfield/
#
#----------------------------------------------------------------------

from skyfield.api import Star, load, Topos
from numpy import arange

class Target():

    def __init__(self):
        # load the planetary positions and our position
        self.planets = load('de421.bsp')
        self.earth = self.planets['earth']
        self.Cas_A = Star(ra_hours=(23, 23, 24.00),
                           dec_degrees=(58, 48, 54.00))

        self.sources = {"Cas_A":self.Cas_A,
                        "Sun":self.planets['sun'],
                        "Moon":self.planets['moon']}

    def get_list(self):
        return list(self.sources)

    def ready_antenna(self, lat, lon, alt):
        self.ant = self.earth + Topos(latitude_degrees=lat,
                        longitude_degrees=lon,
                        elevation_m=alt)

    def ready(self, target):
        self.name = self.sources[target]

    def generate_report(self):
        # the tricky bits to get look angles for a list of times
        ts = load.timescale(builtin=True)

        # get rid of the milliseconds in there
        truth = ts.now().utc
        now_ish = ts.utc(truth[0],truth[1],truth[2],truth[3],truth[4],truth[5])
        now = now_ish.tt

        # make an arry with times incremented by really close to 30 sec
        # this looks weird since it is terrestrial time floating point
        run_time = arange(now, now + .5, .000347222)

        # feed that array to the time function as julian days
        t = ts.tt(jd=run_time)

        # spit out the positon and vector of our location
        astrometric = self.ant.at(t).observe(self.name).apparent()

        # convert this to alt and az lists
        alt, az, distance = astrometric.altaz()

        t_list = t.utc_strftime('%Y:%m:%d:%H:%M:%S')
        alt_list = alt.degrees.tolist()
        az_list = az.degrees.tolist()

        self.report = []
        self.rise_time = None
        self.fade_time = None

        for idx,item in enumerate(t_list):
            if alt_list[idx] < 10.0:
                continue
            else:
                if self.rise_time is None:
                    self.rise_time = t_list[idx]
                self.fade_time = t_list[idx]
                az_out = "{:.4f}".format(az_list[idx])
                alt_out = "{:.4f}".format(alt_list[idx])
                output = t_list[idx]+", "+az_out+", "+alt_out+"\n"
                self.report.append(output)


