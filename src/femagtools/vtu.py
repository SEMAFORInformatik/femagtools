"""Read FEMAG vtu files

"""
import logging
import pathlib
import numpy as np
import vtk
from vtkmodules.util.numpy_support import vtk_to_numpy


class Reader(object):
    '''read vtu-files

        Args:
          pathname : str Directory of result files (vtu-files) or a single vtu file

    '''

    def __init__(self, pathname):
        self.data = {}

        self.reader = vtk.vtkXMLUnstructuredGridReader()
        self.output = self.reader.GetOutput()

        self.field_data_names = []
        self.point_data_names = []
        self.cell_data_names = []
        assert pathlib.Path(pathname).exists(), f"{pathname} not found"
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

        Returns:
            List of values stored in the vtu files

        '''
        return (self.field_data_names +
                self.point_data_names +
                self.cell_data_names)

    def read_data(self, data_list):
        '''Extracts data from the vtu files

        Args:
          data_list : fist of str List of values to extract from vtu_files

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

        Args
          start: float Start of the time window
        end: float End of the time window

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

    def get_field_vector(self, field_data) -> list:
        '''Read field data

        Args:
          field_data : str Name of field to read

        Returns:
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
    def get_point_vector(self, pnt_data, pnt) -> list:
        '''Read point data

        Args:
          point_data : str Name of field to read
          pnt : int Key of point

        Returns:
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

    def get_cell_vector(self, cell_data, cell=0) -> list:
        '''Read cell data

        Args:
          cell_data : str  Name of field to read
          cell : int Key of cell

        Returns:
            List of cell values within the time window

        '''
        if cell_data not in self.data:
            self.read_data([cell_data])

        if cell<=0:
            return vtk_to_numpy(
                self.output.GetCellData().GetAbstractArray(cell_data))
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

    def hrphi(self, scf, elements, magtemp):
        """return H values for each element in polar coord system
        (experimental)
        Args:
          elements: list of model elements
          magtemp: magnet temperature in degree Celsius
          scf: scale factor (poles/poles in model)
        """
        MUE0 = 4e-7*np.pi
        b = vtk_to_numpy(
            self.output.GetCellData().GetAbstractArray("b"))[:, :2]
        e = elements[0]
        x, y = np.mean(vtk_to_numpy(
            self.output.GetCell(e.key-1).GetPoints().GetData()), axis=0)[:2]
        rotpos = np.arctan2(y, x) - np.arctan2(e.center[1], e.center[0])
        psign = 1 if rotpos / scf > 1 else -1

        h = np.empty((len(elements), 2))
        for i, e in enumerate(elements):
            x, y = e.center
            alfa = np.arctan2(y, x)
            b1, b2 = b[e.key-1]
            btempc = 1. + e.br_temp_coef*(magtemp - 20.)
            m1 = e.mag[0]*np.cos(alfa) + e.mag[1]*np.sin(alfa)
            m2 = -e.mag[0]*np.sin(alfa) + e.mag[1]*np.cos(alfa)
            h[i] = abs(e.reluc[0])/MUE0*np.array((
                (b1 - psign*btempc*m1), (b2 - psign*btempc*m2)))
        logging.debug("H shape %s", h.shape)
        return h


    def demag(self, elements):
        """return demag values for each element in kA/m
        Args:
          elements: list of model elements
        Returns:
          list of demags from each vtu file in set
        """
        dlist = []
        for filename in self.filenames:
            self.reader.SetFileName(str(filename))
            self.reader.Update()
            d = vtk_to_numpy(
                self.output.GetCellData().GetAbstractArray("demagnetization"))
            dm = np.empty((len(elements)))
            for i, e in enumerate(elements):
                dm[i] = d[e.key-1]
            dlist.append(dm)
        logging.debug("Demag shape %s", np.array(dlist).shape)
        return dlist


    def get_data_vector(self, data_name, key=0) -> list:
        '''Read data of fiels, point or cell

        Args:
          data_name : str Name of data to read
          key : int (optional) Key of point or cell

        Returns:
            List of values within the time window
        '''
        if data_name in self.field_data_names:
            return self.get_field_vector(data_name)
        if data_name in self.point_data_names:
            return self.get_point_vector(data_name, key)
        if data_name in self.cell_data_names:
            return self.get_cell_vector(data_name, key)
        return []

    def __repr__(self):
        fmt = [f"N Cells {self.output.GetNumberOfCells()}  N Points {self.output.GetNumberOfPoints()}"]
        fmt += [
            f"{self.output.GetPointData().GetAbstractArray(i).GetName()}   Points"
            for i in range(self.output.GetPointData().GetNumberOfArrays())]
        fmt += [
            f"{self.output.GetCellData().GetAbstractArray(i).GetName()}   Cells"
            for i in range(self.output.GetCellData().GetNumberOfArrays())]
        return '\n'.join(fmt)


def read(filename) -> Reader:
    """
    Read vtu file and return Reader object.

    Args:
        filename: name of vtu file to be read
    """
    return Reader(filename)
