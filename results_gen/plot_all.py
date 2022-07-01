# -*- coding: UTF8 -*-
import numpy as np
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime
from netCDF4 import Dataset as Dataset, date2index, num2date
#from grads.gacore import GaCore
# from grads.ganum import GaNum
import argparse
import locale

varname = [None] * 2
levels = [None] * 2
startdate = [None] * 2
enddate = [None] * 2


def parse_args():
    # main =======================================
    parser = argparse.ArgumentParser(description='Evaluate models')

    parser.add_argument("-filename", type=str, help='File name with full path')

    parser.add_argument("-function", type=str,
                        help='function: "bias" for Bias, "avg" for Average, "std" for Standard deviation, "var" for Variance: ')
    parser.add_argument('-plot', dest='isPlot', action='store_true', help='Plot the first case figure on the screen')
    parser.add_argument('-no-plot', dest='isPlot', action='store_false',
                        help='Save the first case figure or function figure, whether specified, to a file')
    parser.set_defaults(isPlot=True)
    parser.add_argument("-savedirectory", type=str, help='Directory to save figures')

    parser.add_argument("-variable", type=str, required=True, help='2D or 3D variable of the first case')
    parser.add_argument("-level", type=str, help='level of the 3D variable  of the first case')
    parser.add_argument("-startdate", type=str, required=True, help='Start date (YYYYMMDDHH)')
    parser.add_argument("-enddate", type=str, required=True, help='End date (YYYYMMDDHH)')

    parser.add_argument("-variable2", type=str, help='2D or 3D variable of the second case')
    parser.add_argument("-level2", type=str, help='level of the 3D variable of the second case')
    parser.add_argument("-startdate2", type=str, help='Start date (YYYYMMDDHH) of the second case')
    parser.add_argument("-enddate2", type=str, help='Forecast date (YYYYMMDDHH) of the second case')

    args = parser.parse_args()
    return args


def opennetcdffile(afilename):
    try:
        print("Opening file ... " + afilename)
        filetoopen = Dataset(afilename, 'r')
    except:
        print("Error opening file ", afilename)
        exit()
    return filetoopen


def plotorsavefigure(isplot, strfile, plt):
    if isplot is False:
        print("Generating figure ... {}".format(strfile))
        plt.savefig(strfile, bbox_inches='tight')
    else:
        print("Viewing figure ... ")
        plt.show()


def drawmaplines(m):
    m.drawcoastlines()
    # m.fillcontinents()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 180., 60.), labels=[0, 0, 0, 1])


def generateFigure(filename, varname, startdate, enddate, savedirectory, isplot, levels=None, y_limits=None, offset=None):

    locale.setlocale(locale.LC_ALL, 'en_US.utf8')

    date_initial = datetime.datetime.strptime(startdate, '%Y%m%d%H')
    date_end = datetime.datetime.strptime(enddate, '%Y%m%d%H')

    # file = opennetcdffile(filename)
    if filename.endswith(".nc"):
        file = opennetcdffile(filename)
        variable = file.variables[varname]
        try:
            varunits = variable.units
        except Exception as e:
            varunits = ' - '

        lat = file.variables['lat'][:]
        lon = file.variables['lon'][:]
        tim = file.variables['time']

        idxdate_ini = date2index(date_initial, tim, calendar='standard')
        idxdate_end = date2index(date_end, tim, calendar='standard')

        dates = num2date(tim[idxdate_ini:idxdate_end + 1], tim.units, calendar='standard')
        if variable.ndim == 3:
            data = variable[idxdate_ini:idxdate_end, 0, 0]
        elif variable.ndim == 4:
            data = variable[idxdate_ini:idxdate_end, :, 0, 0]
            data = np.array(data)
            if levels:
                data = [data[:, i] for i in levels]
                # o loop acima inverte os eixos
                data = np.swapaxes(data, 0, 1)

    plt.figure(figsize=(20, 20), facecolor='white')
    plt.axes().set_xlabel(tim.units)
    if y_limits:
        plt.ylim(y_limits)

    if offset:
        data = np.array(data) + offset

    plt.gca().invert_yaxis()
    # for line, level in zip(data, levels):
    for i in range(data.shape[1]):
        plt.plot(data[:, i], label=levels[i])

    range_x_axis = range(0, data.shape[0], 48)
    plt.xticks(range_x_axis, [dates[di].strftime('%d/%m %h') for di in range_x_axis], rotation=90)
    plt.legend()

    strcaso = "Arquivo: {}, on {}, forecast {} ".format(filename, date_initial.strftime("%d-%m-%Y"), date_end.strftime("%d-%m-%Y"))
    plt.title(strcaso)
    strfile = "{}{}_fore_{}.png".format(savedirectory, filename, enddate)
    plotorsavefigure(isplot, strfile, plt)
    #file.close()

if __name__ == '__main__':

    # -filename /home/denis/_COM_BACKUP/NN_BAM1d/avaliacao/GFCTNMC2014.fct.T062L64___ALL.nc -plot -savedirectory ./
    # -variable temp -startdate 2014011500 -enddate 2014013000 -level 1
    filename = '/home/denis/_COM_BACKUP/NN_BAM1d/avaliacao/IOP2014_CRDRRTMG_HUMO__TIED_GRE_360_savedeach2_64_ORIG.nc'
    savedirectory = "./"
    func = None
    isplot = True


    # variable to avail
    varname[0] = 'temp'
    # levels[0] = list(range(63,0,-1))
    levels[0] = [0, 9, 19, 29, 39, 49, 60, 61, 62, 63]
    limits = [-250, 50]
    offset = -273.15

    # testes
    startdate[0] = '2014010101'
    enddate[0] = '2014033123'

    # # Validating spin up
    # startdate[0] = '2014011501'
    # enddate[0] = '2014013123'

    # Validating umid period
    # startdate[0] = '2014031500'
    # enddate[0] = '2014033123'

    # Validation dry period
    # startdate[0] = '2014091500'
    # enddate[0] = '2014093023'

    if func is None:
        generateFigure(filename, varname[0], startdate[0], enddate[0], savedirectory, isplot, levels=levels[0], y_limits=limits, offset=offset)
    else:
        print("OK")
        #    generateComparison(saveDirectory, isPlot, func, model, mcase, scase, varname, level, dd, mm, yy, hh, ddForecast, mmForecast, yyForecast, hhForecast)


    # filename = args.filename
    # savedirectory = args.savedirectory
    # func = args.function
    # isplot = args.isPlot
    #
    # varname[0] = args.variable
    # level[0] = args.level
    # startdate[0] = args.startdate
    # enddate[0] = args.enddate
    #
    # varname[1] = args.variable2
    # level[1] = args.level2
    # startdate[1] = args.startdate2
    # enddate[1] = args.enddate2

    # This is a simple example which reads a small dummy array, from a
    # netCDF data file created by the companion program simple_xy_wr.py.

    # This example demonstrates the netCDF Python API.
    # It will work either with the Scientific Python NetCDF version 3 interface
    # (http://dirac.cnrs-orleans.fr/ScientificPython/)
    # (http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4_classic-module.html)
    # To switch from one to another, just comment/uncomment the appropriate
    # import statements at the beginning of this file.

    # Jeff Whitaker <jeffrey.s.whitaker@noaa.gov> 20070201

    #  the Scientific Python netCDF 3 interface
    #  http://dirac.cnrs-orleans.fr/ScientificPython/


    # from Scientific.IO.NetCDF import NetCDFFile as Dataset


    # the 'classic' version of the netCDF4 python interface
    # http://code.google.com/p/netcdf4-python/
    # from numpy import arange as # array module from http://numpy.scipy.org
    # from numpy.testing import assert_array_equal, assert_array_almost_equal




    # def generateComparison(saveDirectory, isPlot, strFunction, model, mcase, scase, varname, level, dd, mm, yy, hh, ddForecast, mmForecast, yyForecast, hhForecast):
    #
    #     dateInitial = [None] * 2
    #     dateForecast = [None] * 2
    #     fileName = [None] * 2
    #     ncfile = [None] * 2
    #     variable = [None] * 2
    #     lat = [None] * 2
    #     lon = [None] * 2
    #     tim = [None] * 2
    #     idxDateIni = [None] * 2
    #     idxDateForecast = [None] * 2
    #     dates = [None] * 2
    #     data = [None] * 2
    #     date = [None] * 2
    #
    #     if len(level) != 2:
    #         print("\n\n Error: -level and -level2 paramters are mandatory for comparison! ")
    #         exit()
    #
    #     for i in [0, 1]:
    #
    #         dateInitial[i] = datetime.datetime(yy[i], mm[i], dd[i], hh[i])
    #         dateForecast[i] = datetime.datetime(yyForecast[i], mmForecast[i], ddForecast[i], int(hhForecast[i]))
    #
    #         fileName[i] = getFileName(model[i], mcase[i], scase[i], syy[i], smm[i], sdd[i], shh[i])
    #         ncfile[i] = openNetCdfFile(fileName[i])
    #
    #         try:
    #             variable[i] = ncfile[i].variables[varname[i]]
    #         except Exception:
    #             print("\n\n Error: variable {} not found".format(varname[i]))
    #             exit()
    #
    #         variable[i] = ncfile[i].variables[varname[i]]
    #
    #         lat[i] = ncfile[i].variables['lat'][:]
    #         lon[i] = ncfile[i].variables['lon'][:]
    #         tim[i] = ncfile[i].variables['time']
    #
    #         idxDateIni[i] = date2index(dateInitial[i], tim[i], calendar='standard')
    #         idxDateForecast[i] = date2index(dateForecast[i], tim[i], calendar='standard')
    #         dates[i] = num2date(tim[i][idxDateIni[i]:idxDateForecast[i] + 1], tim[i].units, calendar='standard')
    #         date[i] = dates[i][idxDateForecast[i]]
    #
    #         # in case of one date, use the forecast date ...
    #         if variable[i].ndim == 3:
    #             data[i] = variable[i][idxDateForecast[i], :, :]
    #         elif variable[i].ndim == 4:
    #             data[i] = variable[i][idxDateForecast[i], level[i], :, :]
    #
    #     i = 0
    #     plt.figure(figsize=(20, 20), facecolor='white')
    #     m = Basemap(projection='mill', lat_ts=10, llcrnrlon=lon[i].min(),
    #                 urcrnrlon=lon[i].max(), llcrnrlat=lat[i].min(), urcrnrlat=lat[i].max(),
    #                 resolution='c')
    #
    #     x, y = m(*np.meshgrid(lon[i], lat[i]))
    #
    #     datadiff = np.array(data)
    #     if strFunction == "bias":
    #         strFunctionTitle = "Bias"
    #         datadiff = datadiff[0] - datadiff[1]
    #     elif strFunction == "avg":
    #         strFunctionTitle = "Average"
    #         datadiff = np.mean(datadiff, axis=0)
    #     # elif strFunction == "rmse":
    #         # preciso de mais dados pra calcular as 2 medias
    #         # datadiff = sqrt(pow(media1-media2, 2))
    #     elif strFunction == "std":
    #         strFunctionTitle = "Standard deviation"
    #         datadiff = np.std(datadiff, axis=0)
    #     elif strFunction == "var":
    #         strFunctionTitle = "Variance"
    #         datadiff = np.var(datadiff, axis=0)
    #
    #     m.pcolormesh(x, y, datadiff, shading='flat', cmap=plt.cm.jet)
    #     cbar = m.colorbar(location='bottom', pad="10%")
    #
    #     varUnits = [" - "] * 2
    #     for uni in [0, 1]:
    #         if variable[uni].units != '':
    #             varUnits[uni] = variable[uni].units
    #
    #     levelStr = [None] * 2
    #     if level is None:
    #         levelStr[i] = ""
    #         levelStr[i + 1] = ""
    #     else:
    #         levelStr[i] = ", level = " + str(level[i])
    #         levelStr[i + 1] = ", level = " + str(level[i + 1])
    #
    #     strLabel = varname[i] + r"($" + varUnits[i] + r"$)" + levelStr[i] + " / " + varname[i + 1] + r"($" + varUnits[i + 1] + r"$)" + levelStr[i + 1]
    #     cbar.set_label(strLabel)
    #     drawmaplines(m)
    #
    #     strCaso = strFunctionTitle + ":\nCase 1: {}, {}, {} {}-{}-{} {}:00:00, forecast {}-{}-{} {}:00:00\nCase 2: {}, {}, {} {}-{}-{} {}:00:00, forecast {}-{}-{} {}:00:00" \
    #         .format(model[i], mcase[i], scase[i], syy[i], smm[i], sdd[i], shh[i], syyForecast[i], smmForecast[i], sddForecast[i], shhForecast[i],
    #                 model[i + 1], mcase[i + 1], scase[i + 1], syy[i + 1], smm[i + 1], sdd[i + 1], shh[i + 1], syyForecast[i + 1], smmForecast[i + 1], sddForecast[i + 1],
    #                 shhForecast[i + 1])
    #     plt.title(strCaso)
    #
    #     strFile = "{}{}_{}_{}_{}_{}_{}-{}-{}_{}:00:00_fore_{}-{}-{}_{}:00:00__x__{}_{}_{}_{}_{}-{}-{}_{}:00:00_fore_{}-{}-{}_{}:00:00.png" \
    #         .format(saveDirectory, strFunction, model[i], mcase[i], scase[i], varname[i], syy[i], smm[i], sdd[i], shh[i], syyForecast[i], smmForecast[i], sddForecast[i],
    #                 shhForecast[i], model[i + 1], mcase[i + 1], scase[i + 1], varname[i + 1], syy[i + 1], smm[i + 1], sdd[i + 1], shh[i + 1], syyForecast[i + 1],
    #                 smmForecast[i + 1], sddForecast[i + 1], shhForecast[i + 1])
    #
    #     plotOrSaveFigure(isPlot, strFile, plt)
    #     for i in [0, 1]:
    #         ncfile[i].close()