#!/usr/bin/env python
#
# create a list of points so we can track a target with the new Orbit ACU
#
#   9APR19 - Kyle Eberhart
#   5MAY20 - Kyle Eberhart - extend this to other targets and antennas
#   15DEC20 - Kyle Eberhart - Extend again for TLE targets
#
#   This project made use of Skyfield, http://rhodesmill.org/skyfield/
#
# ----------------------------------------------------------------------

from skyfield.api import load, wgs84, Star, utc
from skyfield import almanac
from numpy import arange, linspace
import numpy as np
import pandas as pd
from datetime import datetime as DT
from datetime import timedelta as TD
from libs import ddor_rsc


class Target():

    def __init__(self, src=None):
        # load the planetary positions and our position
        self.planets = load('de440s.bsp')
        self.ts = load.timescale(builtin=True)
        self.earth = self.planets['earth']

        self.load_radio_src_cat()
        self.Cas_A = Star(ra_hours=(23, 23, 5.00),
                          dec_degrees=(58, 46, 0.00),
                          epoch=1992.5)
        self.sources = {"CAS_A": self.Cas_A,
                        "SUN": self.planets['sun'],
                        "MOON": self.planets['moon'],
                        "JUPITER": self.planets['jupiter barycenter'],
                        "SATURN": self.planets['saturn barycenter'],
                        "VENUS": self.planets['venus barycenter'],
                        "MARS": self.planets['mars barycenter'],
                        "CYG_A": self.cygnus_a,
                        }

        self.default_antenna = (64.97347, 147.50183, 403.6)

        if src is not None:
            self.src = self.sources.get(src, None)
            self.name = self.src

    def ready_antenna(self, lat=None, lon=None, alt=None, ant=None):
        '''Create and use a wgs84 location
        '''
        if lat is not None:
            self.topos_ant = wgs84.latlon(lat, lon, elevation_m=alt)
        if ant is not None:
            self.topos_ant = ant
        else:
            self.topos_ant = wgs84.latlon(self.default_antenna[0],
                                          self.default_antenna[1],
                                          self.default_antenna[2],
                                          )

        self.ant = self.earth + self.topos_ant

    def load_radio_src_cat(self):
        '''Try and load the DSN X-Band radio source catalog
        '''
        with load.open(ddor_rsc.URL) as f:
            self.radio_cat = ddor_rsc.load_dataframe(f)

        self.cygnus_a = Star.from_dataframe(self.radio_cat.loc['3C 405'])

    def visible_targets(self):
        '''A list of the targets visible from this location in 1 hour
        '''
        now_time = DT.utcnow()
        later = TD(hours=1)
        later_time = now_time + later

        tlater = self.ts.from_datetime(later_time.replace(tzinfo=utc))
        # make an array of az el values for plotting
        az_list = []
        el_list = []
        name_list = []
        for name in self.sources:
            alt, az, dist = (self.ant.at(tlater).observe(self.sources[name]).
                             apparent().altaz())
            if alt.degrees < 5.0:
                continue
            az_list.append(az.degrees)
            el_list.append(alt.degrees)
            name_list.append(name)
        points = pd.DataFrame({'name': name_list, 'az': az_list,
                               'el': el_list})
        points['r'] = np.tan((np.pi/4)-(np.radians(points['el'])/2))
        points['theta'] = np.radians(points['az'])
        return points

    def rise_and_fade(self):
        '''A list of the targets visible from this location.
        '''
        now_time = DT.utcnow()
        later = TD(hours=12)
        later_time = now_time + later

        t0 = self.ts.from_datetime(now_time.replace(tzinfo=utc))
        t1 = self.ts.from_datetime(later_time.replace(tzinfo=utc))

        for item in self.sources:
            f = almanac.risings_and_settings(self.planets, self.sources[item],
                                             self.topos_ant)
            t, y = almanac.find_discrete(t0, t1, f)
            if f(t0):
                print(t0.utc_iso(), item, 'is Up')
            for ti, yi in zip(t, y):
                print(ti.utc_iso(), item, 'will Rise' if yi else 'will Set')

    def create_track(self):
        # the tricky bits to get look angles for a list of times

        # get rid of the milliseconds in there
        truth = self.ts.now().utc
        now_ish = self.ts.utc(truth[0], truth[1], truth[2], truth[3], truth[4],
                              truth[5])
        now = now_ish.tt

        # make an arry with times incremented by really close to 30 sec
        # this looks weird since it is terrestrial time floating point
        run_time = arange(now, now + .5, .000347222)

        # feed that array to the time function as julian days
        t = self.ts.tt(jd=run_time)

        # spit out the positon and vector of our location
        self.astrometric = self.ant.at(t).observe(self.name).apparent()

        # convert this to alt and az lists
        alt, az, distance = self.astrometric.altaz()

        self.t_list = t.utc_strftime('%Y:%m:%d:%H:%M:%S')
        self.alt_list = alt.degrees.tolist()
        self.az_list = az.degrees.tolist()

    def generate_report(self):
        '''Create a report in the format used by the OrbitACU to track
        '''
        self.report = []
        self.rise_time = None
        self.fade_time = None

        for idx, item in enumerate(self.t_list):
            if self.alt_list[idx] < 10.0:
                continue
            else:
                if self.rise_time is None:
                    self.rise_time = self.t_list[idx]
                self.fade_time = self.t_list[idx]
                self.az_out = "{:.4f}".format(self.az_list[idx])
                self.alt_out = "{:.4f}".format(self.alt_list[idx])
                output = (self.t_list[idx] + ", " + self.az_out +
                          ", " + self.alt_out + "\n")
                self.report.append(output)

    def add_scan(self, time):
        '''Add a source scan to the already created track, gen a report
        '''



class TLE_Target():

    def __init__(self):
        # load the planetary positions and our position
        self.ts = load.timescale(builtin=True)

    def load_TLE(self, SCID):
        '''Provide a file path to a TLE set to load
        '''
        tle_file = 'data/{}.txt'.format(SCID)
        self.elements = load.tle_file(tle_file)
        self.sat_by_name = {sat.name: sat for sat in self.elements}
        self.sat_by_number = {sat.model.satnum: sat for sat in self.elements}

    def select_sc(self, name=None, number=None):
        '''Select the spacecraft to track, by name or number.
        '''
        if name is not None:
            self.target = self.sat_by_name[name]
        if number is not None:
            self.target = self.sat_by_number[number]
        else:
            self.target = None

    def ready_antenna(self, lat, lon, alt):
        self.ant = wgs84.latlon(lat, lon, elevation_m=alt)

    def list_passes(self):
        '''A list of the passes 1.0 degree over this antennas horizon.
            starting 12 hours prior and ending in 24 hours.
        '''

        self.passes = []

        tnow = self.ts.now()
        # 12 hours prior
        t0 = self.ts.tt(jd=(tnow.tt - .5))
        # 24 hours later
        t1 = self.ts.tt(jd=(tnow.tt + 1))

        t, events = self.target.find_events(self.ant, t0, t1,
                                            altitude_degrees=1.0)

        aos = None
        los = None

        for ti, event in zip(t, events):
            if event == 0:
                aos = ti
            if event == 2:
                los = ti
                if aos is not None:
                    self.passes.append([aos, los])
                aos = None
                los = None

    def pass_angles(self, sc_pass):
        # make an arry with 600 times evenly distributed between AOS and LOS
        run_time = linspace(sc_pass[0].tt, sc_pass[1].tt, 300)

        # feed that array to the time function as julian days
        t = self.ts.tt(jd=run_time)

        # spit out the positon and vector of our location vs the sc
        difference = self.target - self.ant
        topocentric = difference.at(t)

        # convert this to alt and az lists
        alt, az, distance = topocentric.altaz()
        self.angles = np.stack((alt.degrees, az.degrees, run_time), axis=-1)


def create_scan_vectors(magnitude=5, duration=10):
    '''build a vector transform that will cause the antenna to scan the
        source in Az and El, while tracking it across the sky.

        magnitude is the "width" of the scan in degrees
        duration is how log the full scan should take in minutes
    '''
    # 60 seconds per minute, 2 points per second
    points = duration*60*2
    # we will scan in a Rohdonea curve with 5 petals
    k = 5
    alpha = (magnitude/2)
    circle = 2 * np.pi * np.linspace(0, 1, points)
    az = alpha * np.cos(k * circle) * np.cos(circle)
    el = alpha * np.cos(k * circle) * np.sin(circle)

def secant_correct(az, el):
    '''Due to change in elevation the arc length of  azimuth movement changes.
        this function will correct for this change so that the arc length
        is consistent no matter the elevation. For example 1 degree of length
        at the horizon is 1 degree, 1 degree of length at 75 elevation is 3.86
        degrees...

        az is azimuth offset to modify, relative to the horizon.
        el is the elevation to correct for.
    '''
    output = az * (1/np.cos(np.radians(el)))
    return output

if __name__ == "__main__":
    test = TLE_Target()
    test.ready_antenna(64.97347, 147.50183, 403.6)
    test.load_TLE("data/sentinel6a.tle")
#    print(test.sat_by_name)
    test.select_sc(number=46984)
#    print(test.target)
    test.list_passes()
#    print(test.passes)
    test.pass_angles(test.passes[0])
    print(test.angles)
