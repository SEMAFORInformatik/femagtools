"""
Classes for post processing based on vtu-files of created by FEMAG-TS
"""
__author__ = 'werner b. vetter, ronald tanner'

import femagtools.nc
import femagtools.vtu as vtu
import numpy as np
import scipy.integrate as integrate
import warnings
import pathlib
import logging

logger = logging.getLogger('femagtools.asm')


def _read_sections(lines):
    """return list of TS sections

    sections may be surrounded by lines starting with 'Windings/branches'
    Args:
      param lines (list) lines of TS file to read

    Returns:
      list of sections
    """

    section = []
    title = ''
    group = []
    columns = []
    if lines:
        title = lines[0]
        s = 1
        if 'Groups' in lines[1]:
            s = 2
        for line in lines[s:]:
            i = 0
            if 'Windings/branches' in line:
                i = 1
                yield section[i:]
                section = []
            else:
                section.append(line.strip())
    yield section


def read_st(arg, modelname):
    """read TS result files stv, stc, stm
    Arguments:
    arg: (str or Path) name of directory
    modelname: (str) name of model
    """
    r = {}
    if isinstance(arg, str):
        dir = pathlib.Path(arg)
    elif isinstance(arg, pathlib.Path):
        dir = arg
    else:
        logger.error("Invalid arg type %s", type(arg))
        return
    res = dict()
    for ext in ('stv', 'stc', 'stm'):
        lines = (dir / f"{modelname}.{ext}").read_text().split('\n')
        for s in _read_sections(lines):
            if s:
                k = 2
                if ext == 'stm':
                    k = 1
                    labels = [l.lower()
                              for l in s[k-1].split()[::2]][1:]
                elif ext == 'stc':
                    labels = [f'I{l}'
                              for l in s[k-1].split()[1::2]][1:]
                else:
                    labels = [f'U{l}'
                              for l in s[k-1].split()[1::2]][1:]
                m = np.array(
                    [[float(x) for x in l.split()]
                     for l in s[k:] if l]).T
                if 'time' not in res:
                    res['time'] = m[0].tolist()
                for i, k in enumerate(labels):
                    res[k] = m[i+1].tolist()
    return res


def losscoeff_frequency_to_time(B0, f0, c, exp):
    '''Convert Bertotti-coefficient of frequency domain to time domains coefficient
        Parameters
        ----------
        B0 : float
            Base flux density [T]
        f0 : float
            Base freuency [Hz]
        c : float
            Bertotti-coefficient
        exp : float
            Bertotti-exponent
        Return
        -------
        k : float
            Loss coefficient in time domains

        The conversion is only possible for loss-coefficients with
        equal exponent for flux density and frequency,
        as eddy current losses (cw*(B/B0(*)*2*(f/f0)**2) or
        anomalous losses (ce*(B/B0(*)**1.5*(f/f0)**1.5)
        '''
    y, abserr = integrate.quad(lambda x: np.abs(np.cos(2*np.pi*x))**exp, 0, 1)
    return c/(B0**exp*f0**exp)/((2*np.pi)**exp * y)


class TimeRange(object):
    def __init__(self, vtu_data, nc_model):
        '''Read time vector in and generate an equidistant vector if necessary.
           Also the base frequency is determined.
        Parameters
        ----------
        vtu_data : object
            vtu reader
        nc_model: object
        '''

        try:  # FEMAG-TS files
            data_list = ['time [s]']
            vtu_data.read_data(data_list)
            self.vector = vtu_data.get_data_vector('time [s]')
            self.freq = 1/(self.vector[-1]-self.vector[0] +
                           (self.vector[1]-self.vector[0])/2 +
                           (self.vector[-1]-self.vector[-2])/2)
            dt = self.vector[1]-self.vector[0]
            dt_min = 1e32
            self.equidistant = True
            for i in range(len(self.vector)-2):
                dti = self.vector[i+1]-self.vector[i]
                if dt < 0.999*dti or dt > 1.001*dti:
                    self.equidistant = False
                if dti < dt_min:
                    dt_min = dti
            if not self.equidistant:
                numpnt = int((self.vector[-1]-self.vector[0])/dt_min)
                self.vector_equi = np.linspace(self.vector[0],
                                               self.vector[-1],
                                               num=numpnt)
        except:  # FEMAG-DC files
            speed = nc_model.speed
            self.freq = speed/60*nc_model.pole_pairs
            self.equidistant = True


class Losses(object):
    def __init__(self, modelname, dirname):
        '''Loss calculation for FEMAG-TS simulations
        Parameters
        ----------
        dirname : str
            Name of the model (nc-file)
        ncmodel : object

        '''
        self.vtu_data = vtu.read(dirname)
        self.nc_model = femagtools.nc.read(modelname)
        # Read iron losses coefficients
        self.iron_loss_coefficients = self.nc_model.iron_loss_coefficients
        for c in self.iron_loss_coefficients:
            if c['cw_freq_exp'] == c['cw_ind_exp']:
                c['kw'] = losscoeff_frequency_to_time(
                    c['base_induction'],
                    c['base_frequency'],
                    c['cw'], c['cw_freq_exp'])
            else:
                warnings.warn(
                    'Waterfall method not possible, specify parameter kw')
                kw = 0.0

    def ohm_lossenergy_el(self, el, supel):
        '''Ohmic loss energy of an element
        Parameters
        ----------
        el: object
            Element
        supel: object
            Superelement of the element (for material data)

        Returns
        -------
        ellossenergy : float
            Ohmic loss energy of the element
        '''

        length = self.nc_model.arm_length
        time = self.times.vector

        ff = supel.fillfactor
        if ff == 0.0:
            ff = 1.0

        temp_corr = 1+supel.temp_coef*(el.temperature-20)
        ellossenergy = 0.0
        cd_vec = self.vtu_data.get_data_vector('curd', el.key)
        for j in range(len(time)-1):
            cd = (cd_vec[j]+cd_vec[j+1])/2
            dt = time[j+1]-time[j]
            ellossenergy = ellossenergy + dt*cd**2*el.area/ff / \
                supel.conduc*temp_corr*supel.length
        return ellossenergy*length

    def ohm_lossenergy_sr(self, sr):
        '''Ohmic loss energy of a subregion
        Parameters
        ----------
        sr : object
            Subregion

        Returns
        -------
        lossenergy : float
            Ohmic loss energy of the subregion

        The loss energy is determined by adding up the loss energy of the
        individual elements.
        '''
        scale_factor = self.nc_model.scale_factor()
        srlossenergy = 0.0
        for se in sr.superelements:
            selossenergy = 0.0
            if se.conduc > 0.0:
                for el in se.elements:
                    ellossenergy = self.ohm_lossenergy_el(el, se)
                    selossenergy = selossenergy + ellossenergy * scale_factor

            srlossenergy = srlossenergy + selossenergy

        return srlossenergy

    def ohm_lossenergy_subregion(self, srname, start=0.0, end=0.0):
        '''Ohmic loss energy of a subregion
        Parameters
        ----------
        srname:  str
            Name of subregion
        start: float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        lossenergy : float
            Ohmic loss energy of the subregion

        The loss energy is determined by adding up the loss energy of the
        individual elements over the time window.
        If start and end are not specified, the time window of the
        previous calculation is used.
        '''
        data_list = ['time [s]', 'curd']
        self.vtu_data.read_data(data_list)

        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        self.times = TimeRange(self.vtu_data, self.nc_model)

        sr = self.nc_model.get_subregion(srname)
        return self.ohm_lossenergy_sr(sr)

    def ohm_powerlosses_subregion(self, srname, start=0.0, end=0.0):
        '''Ohmic loss dissipation of a subregion within the time window
        Parameters
        ----------
        srname : str
            Name of subregion
        start : float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        powerlosses : float
            Ohmic loss dissipation of the subregion

        The loss energy is determined by adding up the loss energy of the
        individual elements over the time window.
        The loss energy is divided by the time window length
        to obtain the averaged power loss
        If start and end are not specified, the time window of the
        previous calculation is used.
        '''
        while len(srname) < 4:
            srname = srname+' '

        time = self.times.vector[-1]-self.times.vector[0]
        return srlossenergy / time

    def ohm_lossenergy(self, start=0.0, end=0.0):
        '''Ohmic loss energy of all subregions
        Parameters
        ----------
        start: float
            Start of the time window (optional)
        end: float
            End of the time window (optional)

        Returns
        -------
        loss_data: dict
            Dictonary of subregions and ohmic loss energy of it

        The loss energy is determined by adding up the loss energy of the
        individual elements over the time window.
        If start and end are not specified, the time window of the
        previous calculation is used.
        '''

        data_list = ['time [s]', 'curd']
        self.vtu_data.read_data(data_list)

        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        self.times = TimeRange(self.vtu_data, self.nc_model)

        loss_data = []
        for sr in self.nc_model.subregions:
            srlossenergy = self.ohm_lossenergy_sr(sr)

            srname = sr.name
            if sr.wb_key >= 0:
                if srname == '    ':
                    srname = "wdg "+str(sr.wb_key+1)

            loss_data.append(
                {'key': sr.key, 'name': srname, 'losses': srlossenergy})

        return loss_data

    def ohm_powerlosses(self, start=0.0, end=0.0):
        '''Ohmic loss dissipation of all subregions
        Parameters
        ----------
        start : float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        loss_data : dict
            Dictonary of subregions and ohmic loss dissipation of it

        The loss energy is determined by adding up the loss energy of the
        individual elements over the time window.
        The loss energy is divided by the time window length
        to obtain the averaged power loss
        If start and end are not specified, the time window of the
        previous calculation is used.
        '''

        data_list = ['time [s]', 'curd']
        self.vtu_data.read_data(data_list)

        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        self.times = TimeRange(self.vtu_data, self.nc_model)
        #self.time_vector = self.vtu_data.get_data_vector('time [s]')
        time = self.times.vector[-1]-self.times.vector[0]

        loss_data = []
        for sr in self.nc_model.subregions:
            srlossenergy = self.ohm_lossenergy_sr(sr)
            srpowerlosses = srlossenergy / time

            srname = sr.name
            if sr.wb_key >= 0:
                if srname == '    ':
                    srname = "wdg "+str(sr.wb_key+1)

            loss_data.append(
                {'key': sr.key, 'name': srname, 'losses': srpowerlosses})

        return loss_data

    def ohm_powerlosses_fft_el(self, el, supel):
        '''Power dissipation of an element
        Parameters
        ----------
        el : object
            Element
        supel: object
            Superelement of the element (for material data)

        Returns
        -------
        elpowerlosses : float
            Ohmic power losses of the element

        A FFT from the current density is made.
        The power losses of each harmonic is determined and added.
        '''

        length = self.nc_model.arm_length

        ff = supel.fillfactor
        if ff == 0.0:
            ff = 1.0

        temp_corr = 1+supel.temp_coef*(el.temperature-20)
        elpowerlosses = 0.0
        cd_vec_0 = self.vtu_data.get_data_vector('curd', el.key)
        if not self.times.equidistant:
            cd_vec = np.interp(self.times.vector_equi,
                               self.times.vector, cd_vec_0,
                               period=1.0/self.times.freq)
        else:
            cd_vec = cd_vec_0
        cd_spec = abs(np.fft.fft(cd_vec))/(len(cd_vec)/2)
        for j in range(int(len(cd_vec)/2)):
            elpowerlosses = elpowerlosses + \
                cd_spec[j]**2/2*el.area/ff / \
                supel.conduc*temp_corr*supel.length

        return elpowerlosses*length

    def ohm_powerlosses_fft_sr(self, sr):
        '''Power dissipation of a subregion
        Parameters
        ----------
        sr : object
            Subregion

        Returns
        -------
        powerlosses : float
            Ohmic power losses of the subregion

        A FFT from the current density is made.
        The power losses of each harmonic is determined and added.
        '''
        scale_factor = self.nc_model.scale_factor()

        srpowerlosses = 0.0
        for se in sr.superelements:
            sepowerlosses = 0.0
            if se.conduc > 0.0:
                for el in se.elements:
                    elpowerlosses = self.ohm_powerlosses_fft_el(el, se)
                    sepowerlosses = sepowerlosses + elpowerlosses * scale_factor

            srpowerlosses = srpowerlosses + sepowerlosses

        return srpowerlosses

    def ohm_powerlosses_fft_subregion(self, srname, start=0.0, end=0.0):
        '''Power dissipation of a subregion
        Parameters
        ----------
        srname:  str
            Name of subregion
        start: float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        lossenergy : float
            Power dissipation of the subregion

        A FFT from the current density is made.
        The power losses of each harmonic is determined and added.
        The time window has to be pariode or a multiple of it.
        If start and end are not specified, the time window of the
        previous calculation is used.
        '''
        data_list = ['time [s]', 'curd']
        self.vtu_data.read_data(data_list)

        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        self.times = TimeRange(self.vtu_data, self.nc_model)

        sr = self.nc_model.get_subregion(srname)
        srpowerlosses = self.ohm_powerlosses_fft_sr(sr)
        return srpowerlosses

    def ohm_powerlosses_fft(self, start=0.0, end=0.0):
        '''Power dissipation of all subregions
        Parameters
        ----------
        start : float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        loss_data : dict
            Dictonary of subregions and power dissipation of it

        A FFT from the current density is made.
        The power losses of each harmonic is determined and added.
        The time window has to be pariode or a multiple of it.
        If start and end are not specified, the time window of the
        previous calculation is used.
        '''

        data_list = ['time [s]', 'curd']
        self.vtu_data.read_data(data_list)

        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        self.times = TimeRange(self.vtu_data, self.nc_model)

        loss_data = []
        for sr in self.nc_model.subregions:
            srpowerlosses = self.ohm_powerlosses_fft_sr(sr)

            srname = sr.name
            if sr.wb_key >= 0:
                if srname == '    ':
                    srname = "wdg "+str(sr.wb_key+1)

            loss_data.append(
                {'key': sr.key, 'name': srname, 'losses': srpowerlosses})

        return loss_data

# iron losses
    def iron_losses_fft_el(self, el, se):
        '''Iron losses of an element
        Parameters
        ----------
        el: object
            Element
        se: object
            Superelement of element (for material data)

        Returns
        -------
        ironlosses : float
            Iron losses of the element

        A FFT is made from the flux density.
        The iron losses of each harmonic is determined  by
        Bertotti formula

            Physt = ch * (f/f0)**hfe * (B/B0)**hBe * V * rho * shape_factor
            Peddy = ch * (f/f0)**wfe * (B/B0)**wBe * V * rho * shape_factor
            Pexce = ch * (f/f0)**efe * (B/B0)**eBe * V * rho * shape_factor

        and added to the total losses of the element

            Ptot  = (Physt + Peddy + Pexce)
        '''

        length = self.nc_model.arm_length
        freq = self.times.freq

        eltotallosses = 0.0
        elhystlosses = 0.0
        eleddylosses = 0.0
        elexcelosses = 0.0

        if (se.elements[0].reluc[0] < 1.0 or se.elements[0].reluc[1] < 1.0) and \
           (se.elements[0].mag[0] == 0.0 and se.elements[0].mag[1] == 0.0):

            center_pnt = se.elements[0].center
            try:
                if (np.sqrt(center_pnt[0]**2+center_pnt[1]**2) > self.nc_model.FC_RADIUS):
                    ldi = len(self.iron_loss_coefficients)-2  # outside
                else:
                    ldi = len(self.iron_loss_coefficients)-1  # inside
            except:
                ldi = len(self.iron_loss_coefficients)-1  # inside
            #sf = self.iron_loss_coefficients[ldi]['shapefactor']
            if (se.mcvtype > 0):
                ldi = se.mcvtype-1
            bf = self.iron_loss_coefficients[ldi]['base_frequency']
            bb = self.iron_loss_coefficients[ldi]['base_induction']
            ch = self.iron_loss_coefficients[ldi]['ch']
            chfe = self.iron_loss_coefficients[ldi]['ch_freq_exp']
            chbe = self.iron_loss_coefficients[ldi]['ch_ind_exp']
            cw = self.iron_loss_coefficients[ldi]['cw']
            cwfe = self.iron_loss_coefficients[ldi]['cw_freq_exp']
            cwbe = self.iron_loss_coefficients[ldi]['cw_ind_exp']
            ce = self.iron_loss_coefficients[ldi]['ce']
            cefe = self.iron_loss_coefficients[ldi]['ce_freq_exp']
            cebe = self.iron_loss_coefficients[ldi]['ce_ind_exp']
            sw = self.iron_loss_coefficients[ldi]['spec_weight']*1000
            ff = self.iron_loss_coefficients[ldi]['fillfactor']

            bx_vec_0 = self.vtu_data.get_data_vector('b', el.key)[0]
            if not self.times.equidistant:
                bx_vec = np.interp(self.times.vector_equi,
                                   self.times.vector, bx_vec_0,
                                   period=1.0/self.times.freq)
                # f = interpolate.interp1d(self.times.vector, bx_vec_0, kind="cubic")
                # bx_vec = f(self.times.vector_equi)
            else:
                bx_vec = bx_vec_0
            bx_vec = np.array(bx_vec)/ff
            spx = np.fft.fft(bx_vec)
            bx_spec = abs(spx)/(len(bx_vec)/2)
            bx_spec[0] = bx_spec[0]/2
            bx_phi = np.arctan2(spx.imag, spx.real)

            by_vec_0 = self.vtu_data.get_data_vector('b', el.key)[1]
            if not self.times.equidistant:
                by_vec = np.interp(self.times.vector_equi,
                                   self.times.vector, by_vec_0,
                                   period=1.0/self.times.freq)
                # f = interpolate.interp1d(self.times.vector, by_vec_0, kind="cubic")
                # by_vec = f(self.times.vector_equi)
            else:
                by_vec = by_vec_0
            by_vec = np.array(by_vec)/ff
            spy = np.fft.fft(by_vec)
            by_spec = abs(spy)/(len(by_vec)/2)
            by_spec[0] = by_spec[0]/2
            by_phi = np.arctan2(spy.imag, spy.real)

            b_abs = np.sqrt(np.array(bx_vec)**2+np.array(by_vec)**2)
            b_spec = np.sqrt(bx_spec**2+by_spec**2)

            # kh: Korrekturfaktor für drehendes Feld
            kz = 1.0
            i_max = list(b_spec).index(max(b_spec))
            # Phasenverschiebung zwischen Grundschwingungskomponeneten
            dphi = bx_phi[i_max]-by_phi[i_max]
            while dphi > np.pi/2:
                dphi = dphi - np.pi
            while dphi < -np.pi/2:
                dphi = dphi + np.pi
            # Gleichanteil
            b_dc = max(abs(bx_spec[0]), abs(by_spec[0]))
            # Achsenverhältnis der Ellipse
            axis = 0.0
            #Transformation in Hauptrichtung
            j_max = list(b_abs).index(max(b_abs))
            phi = np.arctan2(by_vec[j_max], bx_vec[j_max])
            bxt_vec = []
            byt_vec = []
            for i in range(len(bx_vec)):
                bxt_vec.append(np.cos(phi)*bx_vec[i]+np.sin(phi)*by_vec[i])
                byt_vec.append(np.sin(phi)*bx_vec[i]-np.cos(phi)*by_vec[i])
            # Achsenverhältnis
            max_bxt = max(abs(np.array(bxt_vec)))
            max_byt = max(abs(np.array(byt_vec)))
            if (max_byt > 1.0e-3):
                axis = max_bxt/max_byt
            if axis > 1.0:
                axis = 1.0/axis

            # Korrekturfaktor
            if ((abs(dphi) > np.pi/3) and axis > 0.3):
                kz = 1.55
            if (b_dc > 0.2):
                kz = 1.0 + 0.65*b_dc**2.1
            if (max(b_spec) > 1.85):
                kz = 1.1

            for j in range(int(len(b_spec)/2)):
                elhystlosses = elhystlosses + kz * ch * \
                    (j*freq/bf)**chfe*(b_spec[j]/bb)**chbe
                eleddylosses = eleddylosses + cw * \
                    (j*freq/bf)**cwfe*(b_spec[j]/bb)**cwbe
                elexcelosses = elexcelosses + ce * \
                    (j*freq/bf)**cefe*(b_spec[j]/bb)**cebe

            elhystlosses = elhystlosses*el.area*length*ff*sw
            eleddylosses = eleddylosses*el.area*length*ff*sw
            elexcelosses = elexcelosses*el.area*length*ff*sw

        eltotallosses = elhystlosses + eleddylosses + elexcelosses

        return {'total': eltotallosses,
                'hysteresis': elhystlosses,
                'eddycurrent': eleddylosses,
                'excess': elexcelosses}

    def iron_losses_fft_se(self, se):
        '''Iron losses of a superelement
        Parameters
        ----------
        se: object
            Superelement

        Returns
        -------
        ironlosses : float
            Iron losses of the superelement

        The iron losses are calculated based on the Bertotti formula
        (see also iron_losses_fft_el)
        The results are muliplied by the scale factor,
        represent also the whole machine.
        '''

        scale_factor = self.nc_model.scale_factor()

        sehystlosses = 0.0
        seeddylosses = 0.0
        seexcelosses = 0.0
        if (se.elements[0].reluc[0] < 1.0 or se.elements[0].reluc[1] < 1.0) and \
                (se.elements[0].mag[0] == 0.0 and se.elements[0].mag[1] == 0.0):
            for el in se.elements:
                ellosses = self.iron_losses_fft_el(el, se)
                sehystlosses = sehystlosses + \
                    ellosses['hysteresis'] * scale_factor
                seeddylosses = seeddylosses + \
                    ellosses['eddycurrent'] * scale_factor
                seexcelosses = seexcelosses + ellosses['excess'] * scale_factor

        setotallosses = sehystlosses + seeddylosses + seexcelosses

        return {'total': setotallosses,
                'hysteresis': sehystlosses,
                'eddycurrent': seeddylosses,
                'excess': seexcelosses}

    def iron_losses_fft_subregion(self, srname, start=0.0, end=0.0):
        '''Iron losses of a subregion
        Parameters
        ----------
        srname:  str
            Name of subregion
        start: float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        losses : dict
            Iron losses of the subregion

        The iron losses are calculated based on the Bertotti formula
        (see also ron_losses_fft_se)
        '''
        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        data_list = ['b']
        self.vtu_data.read_data(data_list)
        self.times = TimeRange(self.vtu_data, self.nc_model)

        srtotallosses = 0.0
        srhystlosses = 0.0
        sreddylosses = 0.0
        srexcelosses = 0.0
        sr = self.nc_model.get_subregion(srname)
        for se in sr.superelements:
            selosses = self.iron_losses_fft_se(se)
            srtotallosses = srtotallosses + selosses['total']
            srhystlosses = srhystlosses + selosses['hysteresis']
            sreddylosses = sreddylosses + selosses['eddycurrent']
            srexcelosses = srexcelosses + selosses['excess']

        srlosses = {'subregion': srname,
                    'total': srtotallosses,
                    'hysteresis': srhystlosses,
                    'eddycurrent': sreddylosses,
                    'excess': srexcelosses
                    }
        return srlosses

    def iron_losses_fft(self, start=0.0, end=0.0):
        '''Iron losses of all subregion and superelements
        Parameters
        ----------
        start: float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        losses : dict
            Iron losses of the subregion

        The iron losses are calculated based on the Bertotti formula
        (see also iron_losses_fft_se)
        '''

        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        data_list = ['b']
        self.vtu_data.read_data(data_list)
        self.times = TimeRange(self.vtu_data, self.nc_model)

        losseslist = []

        for se in self.nc_model.superelements:
            selosses = self.iron_losses_fft_se(se)

            if se.subregion:
                for sr in self.nc_model.subregions:
                    if se in sr.superelements:
                        srname = sr.name
            else:
                if (se.mcvtype == 0):
                    center_pnt = se.elements[0].center
                    if (np.sqrt(center_pnt[0]**2+center_pnt[1]**2) > self.nc_model.FC_RADIUS):
                        srname = "no, outside"
                    else:
                        srname = "no, inside"

            found = False
            for srlosses in losseslist:
                if srlosses['subregion'] == srname:
                    srlosses['total'] = srlosses['total']+selosses['total']
                    srlosses['hysteresis'] = srlosses['hysteresis'] + \
                        selosses['hysteresis']
                    srlosses['eddycurrent'] = srlosses['eddycurrent'] + \
                        selosses['eddycurrent']
                    srlosses['excess'] = srlosses['excess']+selosses['excess']
                    found = True
            if not found:
                if selosses['total'] > 0.0:
                    srlosses = {'subregion': srname,
                                'total': selosses['total'],
                                'hysteresis': selosses['hysteresis'],
                                'eddycurrent': selosses['eddycurrent'],
                                'excess': selosses['excess']
                                }
                    losseslist.append(srlosses)

        return losseslist

    def iron_lossenergy_time_el(self, el, se):
        '''Iron losses of an elemt in time domain
        Parameters
        ----------
        el: object
             Element
        se: object
             Superelement of element

        Returns
        -------
        lossenergies : float
             Iron losses of the element

        The iron losses are calculated based on the Bertotti formula
        in time domaine.
        The loss coefficients in frequency domain are converted into
        time domain coefficients.
        For the hysteresis losses is a water fall methode implemented.
        Eddy current losses and anomalous losses are calculated by
        add up the losses of each time step.
        '''

        length = self.nc_model.arm_length
        time = self.times.vector

        elhystenergy = 0.0
        eleddyenergy = 0.0
        elexceenergy = 0.0

        if (se.elements[0].reluc[0] < 1.0 or se.elements[0].reluc[1] < 1.0) and \
                (se.elements[0].mag[0] == 0.0 and se.elements[0].mag[1] == 0.0):
            if (se.mcvtype == 0):
                center_pnt = se.elements[0].center
                try:
                    if (np.sqrt(center_pnt[0]**2+center_pnt[1]**2) > self.nc_model.FC_RADIUS):
                        ldi = len(self.iron_loss_coefficients)-2  # outside
                    else:
                        ldi = len(self.iron_loss_coefficients)-1  # inside
                except:
                    ldi = len(self.iron_loss_coefficients)-1  # use inside
            else:
                ldi = se.mcvtype-1
            kh = self.iron_loss_coefficients[ldi]['kh']
            chbe = self.iron_loss_coefficients[ldi]['ch_ind_exp']
            khml = self.iron_loss_coefficients[ldi]['khml']
            kw = self.iron_loss_coefficients[ldi]['kw']
            cwbe = self.iron_loss_coefficients[ldi]['cw_ind_exp']
            ke = self.iron_loss_coefficients[ldi]['ke']
            cebe = self.iron_loss_coefficients[ldi]['ce_ind_exp']
            sw = self.iron_loss_coefficients[ldi]['spec_weight']*1000
            ff = self.iron_loss_coefficients[ldi]['fillfactor']
            #sf = self.iron_loss_coefficients[ldi]['shapefactor']

            bx_vec = np.array(self.vtu_data.get_data_vector('b', el.key)[0])/ff
            by_vec = np.array(self.vtu_data.get_data_vector('b', el.key)[1])/ff

            # Maximalwert und Richtung des Haupfeldes
            Bpeak = np.sqrt(bx_vec[0]**2+by_vec[0]**2)
            phi = np.arctan2(by_vec[0], bx_vec[0])
            for i in range(1, len(time)):
                b1 = np.sqrt(bx_vec[i-1]**2+by_vec[i-1]**2)
                b2 = np.sqrt(bx_vec[i]**2+by_vec[i]**2)
                if abs(b2) > Bpeak:
                    Bpeak = abs(b2)
                    phi = np.arctan2(by_vec[i], bx_vec[i])

            #Transformation in Hauptrichutng
            br_vec = []
            bt_vec = []
            for i in range(len(time)):
                br_vec.append(np.cos(phi)*bx_vec[i]+np.sin(phi)*by_vec[i])
                bt_vec.append(np.sin(phi)*bx_vec[i]-np.cos(phi)*by_vec[i])

            Bpeak_p = np.sqrt(bx_vec[0]**2+by_vec[0]**2)
            Bx = []
            tp_beg = 0.0
            tp_end = 0.0
            Tp = 0.0
            nzeros = 0
            zero = (br_vec[0] >= 0)
            if br_vec[1] > br_vec[0]:
                up = True
            else:
                up = False
            for i in range(1, len(time)):
                b1 = np.sqrt(br_vec[i-1]**2+bt_vec[i-1]**2)
                b2 = np.sqrt(br_vec[i]**2+bt_vec[i]**2)
                # Maximalwert innerhalb letzter Periode
                if abs(b2) > Bpeak_p:
                    Bpeak_p = abs(b2)
                # Nulldurchgaenge und Periodendauer
                if zero != (br_vec[i] >= 0):
                    zero = (not zero)
                    tp_beg = tp_end
                    tp_end = time[i]
                    if tp_beg > 0.0:
                        nzeros = nzeros+1
                        if nzeros > 1:
                            #Tp = (Tp*(nzeros-1)/nzeros+2*(tp_end-tp_beg)/nzeros)/2
                            Tp = 2*(tp_end-tp_beg)
                            Bpeak = Bpeak_p
                            elhystenergy = elhystenergy+kh*Bpeak**chbe/2
                            Bpeak_p = 0.0
                        else:
                            Tp = 2.0*(tp_end-tp_beg)
                            Bpeak = Bpeak_p
                            elhystenergy = elhystenergy+kh * \
                                Bpeak**chbe * (tp_end-time[0])/Tp
                            Bpeak_p = 0.0
                    Bx = []
                # Wendepunkte
                if up and b2 < b1:
                    Bx.append(b1)
                if not up and b2 > b1:
                    Bx.append(b1)
                # Steigungsrichtung
                if b2 > b1:
                    up = True
                else:
                    up = False

                try:
                    if b2 > 0 and up and b2 > Bx[-2]:
                        Bm = abs(Bx[-2]+Bx[-1])/2
                        dB = abs(Bx[-2]-Bx[-1])
                        elhystenergy = elhystenergy + \
                            kh*Bm**(chbe-1)*khml*dB/2
                        Bx.remove(Bx[-2])
                        Bx.remove(Bx[-1])
                    if b2 < 0 and not up and b2 < Bx[-2]:
                        elhystenergy = elhystenergy + \
                            kh*Bm**(chbe-1)*khml*dB/2
                        Bx.remove(Bx[-2])
                        Bx.remove(Bx[-1])
                    if b2 > 0 and not up and Bx[-1] > Bx[-2]:
                        elhystenergy = elhystenergy + \
                            kh*Bm**(chbe-1)*khml*dB/2
                        Bx.remove(Bx[-2])
                        Bx.remove(Bx[-1])
                    if b2 < 0 and up and Bx[-1] < Bx[-2]:
                        elhystenergy = elhystenergy + \
                            kh*Bm**(chbe-1)*khml*dB/2
                        Bx.remove(Bx[-2])
                        Bx.remove(Bx[-1])

                except:
                    pass

                dt = time[i]-time[i-1]
                dbr = br_vec[i]-br_vec[i-1]
                dbt = bt_vec[i]-bt_vec[i-1]
                db = np.sqrt(dbr**2+dbt**2)
                eleddyenergy = eleddyenergy + kw*(db/dt)**cwbe * dt
                elexceenergy = elexceenergy + ke*(db/dt)**cebe * dt

            #elhystenergy = elhystenergy+kh*Bpeak**chbe * T/(time[-1]-time[0])
            if nzeros >= 1:
                elhystenergy = elhystenergy+kh * \
                    Bpeak**chbe * (time[-1]-tp_end)/Tp

            elhystenergy = elhystenergy*el.area*length*ff*sw
            eleddyenergy = eleddyenergy*el.area*length*ff*sw
            elexceenergy = elexceenergy*el.area*length*ff*sw

        eltotalenergy = elhystenergy + eleddyenergy + elexceenergy

        return {'total': eltotalenergy,
                'hysteresis': elhystenergy,
                'eddycurrent': eleddyenergy,
                'excess': elexceenergy}

    def iron_lossenergy_time_se(self, se):
        '''Iron losses of a superelement in time domain
        Parameters
        ----------
        se: object
            Superelement

        Returns
        -------
        lossenergies : float
            Iron losses of the superelement

        The iron losses are calculated based on the Bertotti formula
        in time domain (see also iron_lossenergy_time_el)
        The results are muliplied by the scale factor,
        represent also the whole machine.
        '''

        scale_factor = self.nc_model.scale_factor()

        sehystenergy = 0.0
        seeddyenergy = 0.0
        seexceenergy = 0.0
        if (se.elements[0].reluc[0] < 1.0 or se.elements[0].reluc[1] < 1.0) and \
                (se.elements[0].mag[0] == 0.0 and se.elements[0].mag[1] == 0.0):
            for el in se.elements:
                elenergy = self.iron_lossenergy_time_el(el, se)
                sehystenergy = sehystenergy + \
                    elenergy['hysteresis'] * scale_factor
                seeddyenergy = seeddyenergy + \
                    elenergy['eddycurrent'] * scale_factor
                seexceenergy = seexceenergy + elenergy['excess'] * scale_factor

        setotalenergy = sehystenergy + seeddyenergy + seexceenergy

        return {'total': setotalenergy,
                'hysteresis': sehystenergy,
                'eddycurrent': seeddyenergy,
                'excess': seexceenergy}

    def iron_lossenergy_time_subregion(self, srname, start=0.0, end=0.0):
        '''Iron loss energy of a subregion
        Parameters
        ----------
        srname:  str
            Name of subregion
        start: float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        losses : dict
            Iron losses energy of the subregion

        The iron losses are calculated based on the Bertotti formula
        in time domain (see also iron_lossenergy_time_el)
        '''
        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        data_list = ['b']
        self.vtu_data.read_data(data_list)
        self.times = TimeRange(self.vtu_data, self.nc_model)

        srtotalenergy = 0.0
        srhystenergy = 0.0
        sreddyenergy = 0.0
        srexceenergy = 0.0
        sr = self.nc_model.get_subregion(srname)
        for se in sr.superelements:
            seenergy = self.iron_lossenergy_time_se(se)
            srtotalenergy = srtotalenergy + seenergy['total']
            srhystenergy = srhystenergy + seenergy['hysteresis']
            sreddyenergy = sreddyenergy + seenergy['eddycurrent']
            srexceenergy = srexceenergy + seenergy['excess']

        return {'subregion': srname,
                'total': srtotalenergy,
                'hysteresis': srhystenergy,
                'eddycurrent': sreddyenergy,
                'excess': srexceenergy
                }

    def iron_losses_time_subregion(self, srname, start=0.0, end=0.0):
        '''Iron power losses of a subregion
        Parameters
        ----------
        srname:  str
            Name of subregion
        start: float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        losses : dict
            Iron losses energy of the subregion

        The iron losses are calculated based on the Bertotti formula
        in time domain (see also iron_lossenergy_time_el)
        '''
        while len(srname) < 4:
            srname = srname+' '

        srenergy = self.iron_lossenergy_time_subregion(srname, start, end)
        time = self.times.vector[-1]-self.times.vector[0]

        srlosses = {'subregion': srname,
                    'total': srenergy['total']/time,
                    'hysteresis': srenergy['hysteresis']/time,
                    'eddycurrent': srenergy['eddycurrent']/time,
                    'excess': srenergy['excess']/time
                    }

        return srlosses

    def iron_lossenergy_time(self, start=0.0, end=0.0):
        '''Iron losses of all subregion and superelements
        Parameters
        ----------
        start: float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        energies : dict
            Iron losses enegies of the subregion

        The iron losses are calculated based on the Bertotti formula
        in time domain (see also iron_lossenergy_time_se)
        '''

        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        data_list = ['b']
        self.vtu_data.read_data(data_list)
        self.times = TimeRange(self.vtu_data, self.nc_model)

        energylist = []

        for se in self.nc_model.superelements:
            selossenergy = self.iron_lossenergy_time_se(se)

            if se.subregion:
                for sr in self.nc_model.subregions:
                    if se in sr.superelements:
                        srname = sr.name
            else:
                if (se.mcvtype == 0):
                    center_pnt = se.elements[0].center
                    try:
                        if (np.sqrt(center_pnt[0]**2+center_pnt[1]**2) > self.nc_model.FC_RADIUS):
                            srname = "no, outside"
                        else:
                            srname = "no, inside"
                    except:
                        srname = "no, used inside"

            found = False
            for srlosses in energylist:
                if srlosses['subregion'] == srname:
                    srlosses['total'] = srlosses['total']+selossenergy['total']
                    srlosses['hysteresis'] = srlosses['hysteresis'] + \
                        selossenergy['hysteresis']
                    srlosses['eddycurrent'] = srlosses['eddycurrent'] + \
                        selossenergy['eddycurrent']
                    srlosses['excess'] = srlosses['excess'] + \
                        selossenergy['excess']
                    found = True
            if not found:
                if selossenergy['total'] > 0.0:
                    srlosses = {'subregion': srname,
                                'total': selossenergy['total'],
                                'hysteresis': selossenergy['hysteresis'],
                                'eddycurrent': selossenergy['eddycurrent'],
                                'excess': selossenergy['excess']
                                }
                    energylist.append(srlosses)

        return energylist

    def iron_losses_time(self, start=0.0, end=0.0):
        '''Iron losses of all subregion and superelements
        Parameters
        ----------
        start: float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------
        losses : dict
            Iron losses of the subregion

        The iron losses are calculated based on the Bertotti formula
        in time domain (see also iron_lossenergy_time_el)
        '''

        energylist = self.iron_lossenergy_time(start, end)
        time = self.times.vector[-1]-self.times.vector[0]

        losseslist = []
        for sr in energylist:
            sr['total'] = sr['total']/time
            sr['hysteresis'] = sr['hysteresis']/time
            sr['eddycurrent'] = sr['eddycurrent']/time
            sr['excess'] = sr['excess']/time
            losseslist.append(sr)

        return losseslist

    def export_lossdensity(self,  filename, methode="fft", start=0.0, end=0.0):
        '''Export the loss density of elements in a vtu -file
        Parameters
        ----------
        filename: string
            Filename of created vtu-file (with extension)
        nethode: string
            Calculation methode (optional, default="fft")
            methode="fft": use fft to calculate the losses
            methode="time": calculate the losses in time domain
        start: float
            Start of the time window (optional)
        end : float
            End of the time window (optional)

        Returns
        -------

        The losses density in each element is calculated und
        stored in a vtu-file.
        '''

        import vtk

        if start != 0.0 or end != 0.0:
            self.vtu_data.set_time_window(start, end)

        #scale_factor = self.nc_model.scale_factor()
        length = self.nc_model.arm_length
        self.times = TimeRange(self.vtu_data, self.nc_model)
        time = self.times.vector[-1]-self.times.vector[0]

        dest_grid = vtk.vtkUnstructuredGrid()
        # copy points
        num_point = self.vtu_data.output.GetNumberOfPoints()
        points = vtk.vtkPoints()
        for i in range(num_point):
            pnt = self.vtu_data.output.GetPoints().GetPoint(i)
            points.InsertNextPoint(pnt)
        dest_grid.SetPoints(points)
        # copy cells
        num_cells = self.vtu_data.output.GetNumberOfCells()
        cells = vtk.vtkCellArray()
        cell_types = []
        # original cells

        for i in range(num_cells):
            cell = self.vtu_data.output.GetCell(i)
            if cell.GetNumberOfPoints() == 3:
                triangle = vtk.vtkTriangle()
                for j in range(cell.GetNumberOfPoints()):
                    triangle.GetPointIds().SetId(j, cell.GetPointId(j))
                cells.InsertNextCell(triangle)
                cell_types.append(vtk.VTK_TRIANGLE)
            if cell.GetNumberOfPoints() == 4:
                quad = vtk.vtkQuad()
                for j in range(cell.GetNumberOfPoints()):
                    quad.GetPointIds().SetId(j, cell.GetPointId(j))
                cells.InsertNextCell(quad)
                cell_types.append(vtk.VTK_QUAD)
        dest_grid.SetCells(cell_types, cells)

        # insert cell values
        cell_data = vtk.vtkDoubleArray()
        cell_data.SetNumberOfComponents(1)
        cell_data.SetName("lossdensity [W/m3]")
        cell_data.SetNumberOfValues(num_cells)

        for se in self.nc_model.superelements:
            for el in se.elements:
                ellossdensity = 0
                if se.conduc > 0.0:
                    if methode == "time":
                        ellossenergy = self.ohm_lossenergy_el(el, se)
                        ellossdensity = ellossenergy / \
                            (time * el.area * length)
                    else:
                        ellosses = self.ohm_powerlosses_fft_el(el, se)
                        ellossdensity = ellosses / (el.area * length)

                if (se.elements[0].reluc[0] < 1.0 or se.elements[0].reluc[1] < 1.0) and \
                   (se.elements[0].mag[0] == 0.0 and se.elements[0].mag[1] == 0.0):
                    if methode == "time":
                        ellossenergy = self.iron_lossenergy_time_el(el, se)
                        ellossdensity = ellossdensity + \
                            ellossenergy['total'] / (time * el.area * length)
                    else:
                        ellosses = self.iron_losses_fft_el(el, se)
                        ellossdensity = ellossdensity + \
                            ellosses['total'] / (el.area * length)

                cell_data.InsertValue(el.key-1, ellossdensity)

        dest_grid.GetCellData().AddArray(cell_data)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(dest_grid)
        writer.Update()
        writer.Write()

        return
