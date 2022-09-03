#!/usr/bin/env python
#
#   Read an antenna.adb file to extract the location and mask info
#
#   6DEC20 - Kyle Eberhart - Created
#
#

#import pandas as pd
import numpy as np
#from scipy.spatial.transform import Rotation as R
import libs.target as T
import libs.link_eng as LE
from skyfield.api import wgs84

TLE_FILE = 'data/sentinel6a.tle'


class Tracker():
    '''Use the antenna class to find information about a target.
    '''

    def __init__(self, file_name, tilt, scid):
        self.ant = Antenna(file_name)
        self.ant.calc_tilt_mask(tilt_az=tilt)
        self.zenith_mask = self.ant.zenith_mask
        self.tilt_mask = self.ant.tilt_mask
        self.mask = self.ant.mask

        self.target = T.TLE_Target()
        self.target.ready_antenna(self.ant.lat, self.ant.lon, self.ant.alt)
        # TODO: make the TLE load more safe and auto
        self.target.load_TLE(scid)
        self.target.select_sc(number=scid)
        self.target.list_passes()
        self.pass_dict = dict(zip([i[0].utc_strftime('%Y/%m/%d-%H:%M:%S')for i
                                   in self.target.passes], self.target.passes))
#        self.pass_angles = pd.DataFrame(columns=['el', 'az', 'color'])

    def calc_pass(self, contact):
        self.target.pass_angles(self.pass_dict[contact])
#        self.pass_angles = pd.DataFrame(self.target.angles,
#                                        columns=['el', 'az', 'color'])


class SunTrack():
    '''Track the sun with an antenna
    '''

    def __init__(self, **kwargs):
        '''
        ant is the the antenna object to use
        start is the start time of the track

        '''
        self.ant = kwargs.get('ant', None)
        self.start = kwargs.get('start', None)
        # sun diam in degrees
        self.sun_diam = .26


class Antenna():
    '''Describe an antenna/ground system
    '''

    def __init__(self, file_name=None, **kwargs):
        '''Read in the file and fill out the variables.

            file_name is the EPOCH mask file, which also has lat, lon and alt.
            diam is antenna diamter in meters
            freq is our targtet frequency in Hz
            lat is latitude in degrees
            lon is longitude in degrees
            alt is altitude in meters
            id is the name of the antenna

        '''
        self.diam = kwargs.get('diam', None)
        self.freq = kwargs.get('freq', None)
        self.lat = kwargs.get('lat', None)
        self.lon = kwargs.get('lon', None)
        self.alt = kwargs.get('alt', None)
        self.id = kwargs.get('id', None)
        if file_name is not None:
            self.file_name = file_name
            self.read_file(self.file_name)
        else:
            self.file_name = None

        if (self.lat is not None) and (self.lon is not None):
            if self.alt is None:
                self.topo = wgs84.latlon(self.lat, self.lon)
            else:
                self.topo = wgs84.latlon(self.lat, self.lon,
                                         elevation_m=self.alt)
            self.loc = self.topo
        if (self.diam is not None) and (self.freq is not None):
            self.set_frequency(self.freq)

    def set_frequency(self, frequency):
        self.freq = frequency
        self.wavelength = LE.calc_wavelength(self.freq)
        # Antennan effiency is a WAG at .67
        self.ant_eff = .528
        self.gain = LE.calc_ant_G(self.ant_eff, self.diam,
                                  self.wavelength)
        self.beamwidth = LE.calc_beamwidth(self.gain)
        self.HPBW = LE.calc_half_power_beamwidth(self.wavelength,
                                                 self.diam)
        self.eff_aperature = LE.calc_effective_aperature(self.gain,
                                                         self.freq)

    def __repr__(self):
        return ("Antenna(file_name={d.file_name!r},id={d.id!r}," +
                "diam={d.diam!r},freq={d.freq!r},lat={d.lat!r}," +
                "lon={d.lon!r},alt={d.alt!r})").format(d=self)

    def __str__(self):
        return self.__repr__()

    def read_file(self, file_name):
        '''You know what it does'''

        self.lat = None
        self.lon = None
        self.alt = None
        self.mask_az = []
        self.mask_el = []
        self.mask = None

        for line in open(file_name):
            if 'position.latitude' in line:
                header, degrees, units = line.split()
                if 'SOUTH' in units:
                    self.lat = float(degrees) * -1
                else:
                    self.lat = float(degrees)
            if 'position.longitude' in line:
                header, degrees, units = line.split()
                if 'WEST' in units:
                    self.lon = float(degrees) * -1
                else:
                    self. lon = float(degrees)
            if 'position.altitude' in line:
                header, altitude, units = line.split()
                if 'KILOMETERS' in units:
                    self.alt = float(altitude) * 1000.0
                else:
                    self.alt = float(altitude)
            if '].azimuth' in line:
                header, azimuth, units = line.split()
                self.mask_az.append(azimuth)
            if '].elevation' in line:
                header, elevation, units = line.split()
                self.mask_el.append(elevation)

#        self.mask = pd.DataFrame({
#            'el': self.mask_el,
#            'az': self.mask_az
#            })

#    def calc_tilt_mask(self, tilt_az=270, offset_angle=7.0):
#        '''Calculate the mask caused by tilting the antenna azimuth off of the
#            horizon so that it can track high angle passes at a lower az
#            velocity.
#        '''
#        # what elevation causes to high of azimuth velocity?
#        zenith_hole = 2.0
#
#        self.tilt_dir = None
#        self.tilt_mask_az = []
#        self.tilt_mask_el = []
#
#        # make a array of vectors for the unit circle in 3d cartesian
#        circle_array = np.linspace(0, 2*np.pi, 359)
#        x = np.cos(circle_array)
#        y = np.sin(circle_array)
#        z = np.zeros(359)
#        vect_array = np.stack((x, y, z), axis=-1)
#
#        # rotate the tilt offset in the x azis
#        vect_tilt = self.rotate(vect_array, np.radians(offset_angle), 'y')
#
#        # rotate the to the tilt angle
#        vect_rot = self.rotate(vect_tilt, np.radians(tilt_az), 'z')
#
#        self.alt_az = self.vect_angle(vect_rot)
#
#        self.tilt_mask = pd.DataFrame(self.alt_az, columns=['el', 'az'])
#
#        # calculate the key hole at zenith, given this tilt info
#        zenith_circle = np.linspace(0, 2*np.pi, 359)
#        # figure the xy for each az angle in zenith_circle
##        x = np.multiply(np.sin(np.radians(zenith_hole)), np.cos(zenith_circle))
#        y = np.multiply(np.sin(np.radians(zenith_hole)), np.sin(zenith_circle))
#        z = np.cos(np.full(359, np.radians(zenith_hole)))
#        zenith_array = np.stack((x, y, z), axis=-1)
#
#        # rotate the tilt offset in the x axis
#        zenith_tilt = self.rotate(zenith_array, np.radians(offset_angle), 'y')
#
#        # rotate to the tilt angle
#        zenith_rot = self.rotate(zenith_tilt, np.radians(tilt_az), 'z')
#
#        self.zenith_alt_az = self.vect_angle(zenith_rot)
#
#        self.zenith_mask = pd.DataFrame(self.zenith_alt_az,
#                                        columns=['el', 'az'])

#    def rotate(self, vect, rot_radians, axis='x'):
#        '''rotate the antenna unit circle around the given axis...
#        '''
#        if axis == 'x':
#            rot_axis = np.array([1, 0, 0])
#        if axis == 'y':
#            rot_axis = np.array([0, 1, 0])
#        if axis == 'z':
#            rot_axis = np.array([0, 0, 1])
#        rot_vector = rot_radians * rot_axis
#        rotation = R.from_rotvec(rot_vector)
#        return rotation.apply(vect)
#
    def vect_angle(self, vect):
        '''Find the angles of the vectors
        '''
        tilt_array = np.empty([359, 2])
        xy = vect[:, 0]**2 + vect[:, 1]**2
        tilt_array[:, 0] = np.rad2deg(np.arctan2(vect[:, 2], np.sqrt(xy)))
        tilt_array[:, 1] = np.rad2deg(np.arctan2(vect[:, 1], vect[:, 0]))
        return tilt_array

    def _fix_angle(self, angle):
        '''Check if an angle is greater than 360 or less than zero
        '''
        if angle > 360:
            return (angle - 360)
        if angle < 0:
            return (360 + angle)
        return angle


if __name__ == '__main__':
    #    test = Antenna('data/fbfe3.antenna.adb')
    #    print(test.mask)
    test = Antenna('data/fbfe3.antenna.adb')
    test.calc_tilt_mask()
    #    print(test.vect)
    print(test.alt_az)
    pass
