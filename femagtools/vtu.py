"""
    femagtools.vtu
    ~~~~~~~~~~~~~~

    Read FEMAG vtu files
"""
import vtk
import pathlib
import numpy as np


class Reader(object):
    '''Class to read vtu-files'''

    def __init__(self, pathname):
        '''Read the vtu-files
        Parameters
        ----------
        pathname : str
            Directory of result files (vtu-files) or a single vtu file
        '''

        self.data = {}

        self.reader = vtk.vtkXMLUnstructuredGridReader()
        self.output = self.reader.GetOutput()

        self.field_data_names = []
        self.point_data_names = []
        self.cell_data_names = []
        if pathlib.Path(pathname).suffix == '.vtu':
            self.filenames = [pathlib.Path(pathname)]
        else:
            self.filenames = sorted(pathlib.Path(pathname).glob(
                "*.vtu"))
        self.reader.SetFileName(str(self.filenames[0]))
        self.reader.Update()
        self.field_data_names = [
            self.output.GetFieldData().GetAbstractArray(i).GetName()
            for i in range(self.output.GetFieldData().GetNumberOfArrays())]
        self.point_data_names = [
            self.output.GetPointData().GetAbstractArray(i).GetName()
            for i in range(self.output.GetPointData().GetNumberOfArrays())]
        self.cell_data_names = [
            self.output.GetCellData().GetAbstractArray(i).GetName()
            for i in range(self.output.GetCellData().GetNumberOfArrays())]

        self.set_time_window(0.0, 0.0)

    def get_data_names(self):
        '''Read the list of values stored in the vtu files
        Parameters
        ----------
            None
        Returns
        -------
            List of values stored in the vtu files
        '''
        return (self.field_data_names +
                self.point_data_names +
                self.cell_data_names)

    def read_data(self, data_list):
        '''Extracts data from the vtu files
        Parameters
        ----------
        data_list : fist of str
            List of values to extract from vtu_files
        '''
        for data_name in data_list:
            if data_name in self.field_data_names:
                self.data[data_name] = []
            elif data_name in self.point_data_names:
                self.data[data_name] = []
            elif data_name in self.cell_data_names:
                self.data[data_name] = []
            else:
                raise Exception('unknown data name "' + data_name+'"')

        for filename in self.filenames:
            self.reader.SetFileName(str(filename))
            self.reader.Update()

            for data_name in data_list:
                if data_name in self.field_data_names:
                    self.data[data_name].append(
                        self.output.GetFieldData().GetAbstractArray(data_name))
                if data_name in self.point_data_names:
                    self.data[data_name].append(
                        self.output.GetPointData().GetAbstractArray(data_name))
                if data_name in self.cell_data_names:
                    self.data[data_name].append(
                        self.output.GetCellData().GetAbstractArray(data_name))

        return "done"

    def set_time_window(self, start, end):
        '''Set time window
        Parameters
        ----------
        start: float
            Start of the time window
        end: float
            End of the time window

        Only values within the time window are output by the functions
            get_field_vector
            get_point_vector
            get_cell_vector
            get_data_vector
        At start = 0.0 the values are read out starting from the first value
        At end = 0.0 the values are read out up to the last value
        '''
        try:
            if "time [s]" not in self.data:
                self.read_data(['time [s]'])

            if start == 0 or start <= self.data['time [s]'][0].GetValue(0):
                self.istart = 0
            else:
                self.istart = 0
                for i in range(len(self.data['time [s]'])):
                    if start >= self.data['time [s]'][i].GetValue(0):
                        self.istart = i+1
            if end == 0 or end >= self.data['time [s]'][-1].GetValue(0):
                self.iend = len(self.data['time [s]'])
            else:
                self.iend = 0
                for i in range(len(self.data['time [s]'])):
                    if end <= self.data['time [s]'][i].GetValue(0):
                        self.iend = i
        except:
            self.istart = None
            self.iend = None

    def get_field_vector(self, field_data):
        '''Read field data
        Parameters
        ----------
        field_data : str
            Name of field to read
        Returns
        -------
        field_vec : list of float   
            List of field values within the time window
        '''
        if field_data not in self.data:
            self.read_data([field_data])

        if self.istart:
            start = self.istart
        else:
            start = 0
        if self.iend:
            end = self.iend
        else:
            end = len(self.data[field_data])

        field_vec = []
        # for i in range(self.istart,self.iend):
        for i in range(start, end):
            field_vec.append(self.data[field_data][i].GetValue(0))
        return field_vec

    # pnt = node-key, >0
    def get_point_vector(self, pnt_data, pnt):
        '''Read point data
        Parameters
        ----------
        point_data : str
            Name of field to read
        pnt : int
            Key of point
        Returns
        -------
        point_vec : list of float   
            List of point values within the time window
        '''
        if pnt_data not in self.data:
            self.read_data([pnt_data])

        if self.istart:
            start = self.istart
        else:
            start = 0
        if self.iend:
            end = self.iend
        else:
            end = len(self.data[pnt_data])

        point_vec = []
        for i in range(start, end):
            point_vec.append(self.data[pnt_data][i].GetValue(pnt-1))
        return point_vec

    def get_cell_vector(self, cell_data, cell):
        '''Read cell data
        Parameters
        ----------
        cell_data : str
            Name of field to read
        cell : int
            Key of cell
        Returns
        -------
        cell_vec : list of float   
            List of cell values within the time window
        '''
        if cell_data not in self.data:
            self.read_data([cell_data])

        i = self.cell_data_names.index(cell_data)
        noc = self.output.GetCellData().GetAbstractArray(i).GetNumberOfComponents()
        if noc == 1:
            cell_vec = []
        else:
            cell_vec_x = []
            cell_vec_y = []
            cell_vec_z = []

        if self.istart:
            start = self.istart
        else:
            start = 0
        if self.iend:
            end = self.iend
        else:
            end = int(len(self.data[cell_data]))

        for i in range(start, end):
            if noc == 1:
                cell_vec.append(self.data[cell_data][i].GetValue(cell-1))
            else:
                cell_vec_x.append(
                    self.data[cell_data][i].GetValue(noc*(cell-1)))
                cell_vec_y.append(
                    self.data[cell_data][i].GetValue(noc*(cell-1)+1))
                cell_vec_z.append(
                    self.data[cell_data][i].GetValue(noc*(cell-1)+2))
        if noc == 1:
            return cell_vec
        else:
            return [cell_vec_x, cell_vec_y, cell_vec_z]

    def get_data_vector(self, data_name, key=0):
        '''Read data of fiels, point or cell
        Parameters
        ----------
        data_name : str
            Name of data to read
        hey : int (optional)
            Key of point or cell
        Returns
        -------
        data_vec : list of float   
            List of values within the time window
        '''
        if data_name in self.field_data_names:
            return self.get_field_vector(data_name)
        if data_name in self.point_data_names:
            return self.get_point_vector(data_name, key)
        if data_name in self.cell_data_names:
            return self.get_cell_vector(data_name, key)
        return []


def read(filename):
    """
    Read vtu file and return Reader object.

    Arguments:
        filename: name of vtu file to be read
    """
    return Reader(filename)
