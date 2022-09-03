#!/usr/bin/env python3
#
#   A library for satellite link engineering equations.
#
#   Kyle Eberhart - 29OCT20
#
#   Public Domain - no implied warranty, use at your own risk
#
#   1. This library uses many of the ideas and equations from;
#       "Link Performance Analysis for a Proposed Future Architecture of the
#       Air Force Satellite Control Network" by Eric W. Nelson, USAF
#
#   2. Just like Cpt. Nelson I also referenced;
#       "DSN Telecommunications Link Design Handbook" 810-005 rev. E
#
#   3. I have also used;
#       "Electromagnetic Waves and Antennas" by Sophocles J. Orfanidis
#
#   4. Noise figure calculations and antenna effective aperature were taken
#       from the source code of "Virgo: A Versatile Spectrometer for Radio
#       Astronomy" by Apostolos Spanakis-Misirlis, Cameron L. Van Eck and
#       E.p. Boven
#
# -----------------------------------------------------------------------------

import math
#from scipy.special import jv
#from scipy.special import erfc

k = 1.3806*math.pow(10, -23)    # J/K
k_dBW = -228.5991               # dBW/K/Hz
c = 299792458                   # m/s
earth_radius = 6378000          # meters


def calc_SNR(EIRP, L, GoT):
    '''Signal to Noise Ratio

        (C/No) = (EIRP)*(1/L)*(GoT)*(1/k)

        EIRP is the Effective Isotropic Radiated Power
        L is the medium losses
        GoT is the receiving system G/T or System Gain over System noise
            temperature.
        k is Boltzmann's constant (1.3806x10^-23 J/K or -228.5991 dBW/K/Hz)

    '''
    _SNR = (EIRP)*(1/L)*(GoT)*(1/k)
    return _SNR


def calc_EIRP(G, P):
    '''Effective Isotropic Radiated Power

        EIRP = G*P

        G is the Gain of the transmit antenna in dB
        P is the radiated power in mW

    '''
    G = math.pow(10, G/10)
    _EIRP = G*P
    return _EIRP


def calc_wavelength(freq):
    '''Wavelength

        wave_length = c/freq

        c is the speed of light in m/s
        freq is the frequency in Hz

    '''
    wavelength = c/freq
    return wavelength


def calc_ant_G(antenna_effiency, diameter, wavelength):
    '''Gain of a simple prime focus parbolic antenna

        G = antenna_effiency*(pi()*diameter/wavelength)^2

        antenna_effiency is a decimal percentage
        diameter is in m
        wavelength is in m

    '''
    _G = antenna_effiency*math.pow(math.pi*diameter/wavelength, 2)
    return 10*math.log(_G, 10)


def calc_effective_aperature(G, freq):
    '''Antenna effective aperature [m^2]

        G is antenna gain in dB
        freq is in Hz

    '''
    A_eff = (10**(G/10)*calc_wavelength(freq)**2)/(4*math.pi)
    return A_eff


def calc_beamwidth(G):
    '''Antenna beamwidth in degrees

        G is antenna gain in dB, calculated above...

    '''
    G = math.pow(10, G/10)
    _beamwidth = math.sqrt(16/G)
    return math.degrees(_beamwidth)


def calc_half_power_beamwidth(wavelength, diameter):
    '''3dB beamwidth or HPBW

        reference 3; page 748; equation 16.3.11
        "The constant 70 degrees represents only a rough approximation..."

        HPBW = 70*(wavelength/diameter)

        diameter of the parabolic reflector in meters
        wavelength is the wavelength in meters

    '''
    _HPBW = 70*wavelength/diameter
    return _HPBW


def calc_antenna_T(beamwidth, antenna_effiency, sky_temp_K, ambient_temp_K):
    '''Antenna Temperature

        beamwidth of the antenna at frequency
        antenna_effiency is a decimal percentage
        sky_temp_K is the sky temperature in Kelvin, varies per frequency
        ambient_temp_K is the ambient temperature in Kelvin

        I am not sure where I got this equation/function. It has parts that
        look similar to things in reference 3, at the end of page 758.

        I would like to update this to better account for frequency and
        elevation angle. Reference 2 has several interesting equations but
        will require quite a bit of work to get things working.

    '''
    Ta_mb = 1/beamwidth*(sky_temp_K*(antenna_effiency)*beamwidth)
    Ta_gbl = 1/beamwidth*(ambient_temp_K*(1-antenna_effiency)/2*beamwidth)
    Ta_hbl = 1/beamwidth*(ambient_temp_K/2*(1-antenna_effiency)/2*beamwidth)
    _T = Ta_mb + Ta_gbl + Ta_hbl
    return _T


def calc_G_T(G, T_sys):
    '''Calculate the antenna G/T

        G is the antenna gain in dBi
        T_sys is the antten noise temperature in K

    '''
    _G_T = G-lin_to_db(T_sys)
    return _G_T


def NF_to_T_noise(NF, T_ref=290):
    '''Convert Noise Figure to noise temperature [K]

        NF is Noise Figure in dB
        T_ref is the reference temperature in K

    '''
    _T_noise = T_ref*((10**(NF/10)) - 1)
    return _T_noise


def T_noise_to_NF(T_noise, T_ref=290):
    '''Convert a noie temperature to NF [dB]

        T_noise is the noise temperature in K
        T_ref is the reference temperature in K

    '''
    _NF = lin_to_db((T_noise/T_ref) +1)
    return _NF


def calc_SEFD(eff_aperature, T_sys):
    '''System Equivilent flux density [Jy]

        eff_aperature is the antenna effective aperature in m^2
        T_sys is the system noise temperature in [K]

    '''
    _sefd = 10**26 * 2*k*T_sys/eff_aperature
    return _sefd


def calc_radiometer_equation(S_flux, sefd, on_time, bw):
    '''Estimate the snr for an observation

        S_flux is the source flux density in Jy
        sefd is the system equivalent flux density in Jy
        on_time is the on source integration time in seconds
        bw is the aquisition bandwidth in Hz

    '''
    _snr = S_Flux*math.sqrt(on_time*bw)/sefd
    return _snr


def calc_off_nadir(el_angle, sc_alt, gs_alt):
    '''Angle off of earth pointng beam

        elevation angle of the ground station antenna
        spacecraft altitude in meters
        ground station altitude in meters

    '''
    gs_part = (gs_alt+earth_radius)*math.sin(math.radians(el_angle+90.0))
    _DOFF = math.degrees(math.asin(gs_part/(sc_alt+earth_radius)))
    return _DOFF


def calc_pointing_loss(point_err, HPBW):
    '''Pointing error loss

        point_err is the pointing offset error in degrees
        HPBW is the half power beamwidth of the antenna in degrees

    '''
    print("doesn't work right!")
    loss1 = 3*math.pow(point_err/HPBW, 2)
    loss2 = -12*math.pow(point_err/HPBW, 2)
    loss3 = 0.063*(math.pow(point_err, 2)/math.pow(HPBW, 2))
    _Point_loss = (10*math.log(
                        math.exp(
                            2.773*math.pow(
                                point_err, 2)/math.pow(HPBW, 2))))
    return _Point_loss, loss1, loss2, loss3


def calc_atmo_loss(el_angle):
    '''A rough estimate for atmospheric loss. There are much better ways
        to do this, but this is very easy.

    '''
    _atmo_loss = 1+(1/el_angle)
    return lin_to_db(_atmo_loss)


def calc_free_space_loss(slant_range, frequency):
    '''Free Space Loss

        FSL = (4*pi*range*frequency/c)^2

        slant_range is the straight line distance to the spacecraft from
            the earth terminal in meters
        frequency is in Hz
        c is the speed of light in m/s

    '''
    _FSL = math.pow(4*math.pi*slant_range*frequency/c, 2)
    return lin_to_db(_FSL)


def calc_polarization_loss(nadir_off):
    '''Polarization Loss

        Lpol = 1.389*10^8(nadir_off^4)-3.389*10^4(nadir_off^2)-2.86*10^7 (dB)

        nadir_off is radians offset from boresight

    '''
    print("doesn't work right!")
    _Lpol = 1.389*math.pow(
            10, -8)*math.pow(
                    nadir_off, 4)-3.389*math.pow(
                        10, -4)*math.pow(
                                nadir_off, 2)-2.286*math.pow(10, -7)
    return _Lpol


def uplink_performance(EIRP, uplink_loss, GoTsc, k):
    '''Uplink performance

        (C/No)uplink = (EIRPground_station)(1/uplink_loss)(G/Tsc)(1/k)

        EIRPground_station is the EIRP from the ground station
        uplink_loss is the total uplink losses in dB
        k is Boltzmann's constant

    '''
    _CoNo = EIRP*(1/uplink_loss)*(GoTsc)*(1/k)
    return _CoNo


def downlink_performance(EIRP, downlink_loss, GoTgs, k):
    '''Downlink performance

        (C/No)downlink = (EIRPsc)*(1/downlink_loss)*(GoTgs)*(1/k)

        EIRPsc is the EIRP from the spacecraft
        downlink_loss is the total downlnk losses in dB
        k is Boltzmann's constant

    '''
    _CoNo = EIRP*(1/downlink_loss)*(GoTgs)*(1/k)
    return _CoNo


#def service_mod_loss(mod_index):
#    '''Service Modulation loss
#
#        service_mod_loss = 10*log10(2*bessel(1, mod_index)^2)
#
#        mod_index is the modulation index of the subcarrier
#        bessel is the bessel function of the first order
#
#    '''
#    _service_mod_loss = 10*math.log(2*math.pow(jv(1, mod_index), 2))
#    return _service_mod_loss


def TLM_EbNo(CoNoTLM, service_mod_loss, data_rate_loss):
    '''Eb/No of the telemetry stream

        CoNoTLM is C/No of the data stream
        service_mod_loss is the loss due to modulation
        data_rate_loss is the losses due to the datarate changes

    '''
    _TLM_EbNo = CoNoTLM-service_mod_loss-data_rate_loss
    return _TLM_EbNo


#def bit_error_rate(TLM_EbNo):
#    '''Bit Error Rate - if it is an SGLS waveform
#
#        BER = 0.5*erfc(sqrt(TLM_EbNo)
#
#        erfc is the complimentary error function
#        TLM_EbNo is the telemetry subcarrier engergy per bit over noise
#            density in dB
#
#    '''
#    _BER = 0.5*erfc(math.sqrt(TLM_EbNo))
#    return _BER
#

def lin_to_db(value):
    '''Convert a value from linear form to logrithmic dB
    '''
    result = 10*math.log10(value)
    return result


def db_to_lin(value):
    '''Convert a value from dB to linear form
    '''
    result = math.pow(10, value/10)
    return result
