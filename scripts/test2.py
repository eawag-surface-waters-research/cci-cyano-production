import os
import numpy as np
import scipy.io as sio
import datetime
from netCDF4 import Dataset


def store_time_series_and_null_info_CCI(lakeNR, lakeMetaData, varName, minDateNR_FULL, maxDateNR_FULL,
                                        replaceTimeSeriesAndNullInfo):
    # Folder to store time series and null info data
    folderNameStoreTimeSeriesAndNullInfoData = f"data/CCISTEP1_timeSeriesAndNullInfo/Lake_{lakeNR}_{lakeMetaData['id'][lakeNR]}_{lakeMetaData['short_name'][lakeNR]}"

    if not os.path.exists(folderNameStoreTimeSeriesAndNullInfoData):
        os.makedirs(folderNameStoreTimeSeriesAndNullInfoData)

    fileName_timeSeriesNullInfo = f"{folderNameStoreTimeSeriesAndNullInfoData}/CCI_NullInfo_{varName}_{minDateNR_FULL}_{maxDateNR_FULL}"
    fileName_timeSeriesData = f"{folderNameStoreTimeSeriesAndNullInfoData}/CCI_Data_{varName}_{minDateNR_FULL}_{maxDateNR_FULL}"

    # Check if file exists and load data if replace is false
    if not replaceTimeSeriesAndNullInfo:
        if os.path.isfile(f"{fileName_timeSeriesNullInfo}.mat"):
            timeSeriesNullInfoData = sio.loadmat(f"{fileName_timeSeriesNullInfo}.mat")
            dataFiles_allFileNames = timeSeriesNullInfoData['dataFiles_allFileNames']
            dataFiles_NRdataFiles = timeSeriesNullInfoData['dataFiles_NRdataFiles']
            print(
                f"storeTimeSeriesAndNullInfo_MerisModis done {lakeMetaData['name'][lakeNR]} ({lakeMetaData['short_name'][lakeNR]}): - Loaded from file")
            return fileName_timeSeriesNullInfo, dataFiles_allFileNames, dataFiles_NRdataFiles

    # Set QA flag name based on variable name
    if varName == 'chla_mean':
        QAflagName = 'lwlr_quality_flag'
    else:
        raise ValueError("could not set QAflag name")

    NRtimePoints = maxDateNR_FULL - minDateNR_FULL + 1
    timeSeries_dateNR = np.nan * np.ones(NRtimePoints)
    timeSeries_DOY = np.nan * np.ones(NRtimePoints)
    timeSeries_dataMatrix = []
    timeSeries_QAMatrix = []
    stats_nrDataPerDOY = np.zeros(366)
    stats_YRs = np.arange(2002, 2023)
    stats_nrDataPerYR = np.zeros(len(stats_YRs))

    timePointNR = 0
    firstDateWithData = True

    for dateNR in range(minDateNR_FULL, maxDateNR_FULL + 1):
        date = datetime.datetime.strptime(str(dateNR), "%Y%m%d")

        # Load chla data
        ncFileName = f"data/CCI_lakes_nc/v2.1.0/{date.year}/{date.strftime('%Y-%m-%d')}"
        matFolder = f"data/CCI_lakes_mat/v2.1.0/{varName}/{date.year}"
        matFileName = f"{matFolder}/ESACCI-LAKES-L3S-LK_PRODUCTS-MERGED-{date.strftime('%Y-%m-%d')}-fv2.1.0.mat"

        if not os.path.exists(matFolder):
            os.makedirs(matFolder)

        timePointNR += 1
        if os.path.isfile(ncFileName):

            # Load data from mat file
            if os.path.isfile(matFileName):
                matFileData = sio.loadmat(matFileName)

                ncVarData_nrRows = matFileData['ncVarData_nrRows']
                ncVarData_nrCols = matFileData['ncVarData_nrCols']
                ncVarData_coorsWithData = matFileData['ncVarData_coorsWithData']
                ncVarData_dataAtCoors = matFileData['ncVarData_dataAtCoors']
                ncQAData_dataAtCoors = matFileData['ncQAData_dataAtCoors']

                ncVarData = np.full((ncVarData_nrRows, ncVarData_nrCols), np.nan)
                ncVarData[ncVarData_coorsWithData] = ncVarData_dataAtCoors
                ncQAData = np.full((ncVarData_nrRows, ncVarData_nrCols), np.nan)
                ncQAData[ncVarData_coorsWithData] = ncQAData_dataAtCoors

                ncLon = matFileData['ncLon']
                ncLat = matFileData['ncLat']

            else:
                # Load from nc file
                ncData = Dataset(ncFileName, 'r')
                ncVarData = ncData.variables[varName][:]
                ncQAData = ncData.variables[QAflagName][:]

                ncVarData_nrRows, ncVarData_nrCols = ncVarData.shape
                ncVarData_coorsWithData = np.where(~np.isnan(ncVarData))
                ncVarData_dataAtCoors = ncVarData[ncVarData_coorsWithData]
                ncQAData_dataAtCoors = ncQAData[ncVarData_coorsWithData]

                ncLon = ncData.variables['lon'][:]
                ncLat = ncData.variables['lat'][:]

                sio.savemat(matFileName, {'ncVarData_nrRows': ncVarData_nrRows, 'ncVarData_nrCols': ncVarData_nrCols,
                                          'ncVarData_coorsWithData': ncVarData_coorsWithData,
                                          'ncVarData_dataAtCoors': ncVarData_dataAtCoors,
                                          'ncQAData_dataAtCoors': ncQAData_dataAtCoors, 'ncLon': ncLon, 'ncLat': ncLat})
                print(f"global {varName} data loaded from nc file and stored as mat file")

            if firstDateWithData:
                lon_min_box = lakeMetaData['lon_min_box'][lakeNR]
                lon_max_box = lakeMetaData['lon_max_box'][lakeNR]
                lat_min_box = lakeMetaData['lat_min_box'][lakeNR]
                lat_max_box = lakeMetaData['lat_max_box'][lakeNR]

                coorsInBox_yDir = np.where((ncLon < lon_max_box) & (ncLon > lon_min_box))
                coorsInBox_xDir = np.where((ncLat < lat_max_box) & (ncLat > lat_min_box))

                firstDateWithData = False

            ncVarData = ncVarData[coorsInBox_yDir, coorsInBox_xDir]
            ncQAData = ncQAData[coorsInBox_yDir, coorsInBox_xDir]

            dateYR = date.year
            DOY = (dateNR - datetime.datetime(dateYR - 1, 12, 31).toordinal())

            timeSeries_dateNR[timePointNR - 1] = dateNR
            timeSeries_DOY[timePointNR - 1] = DOY

            if not timeSeries_dataMatrix:
                timeSeries_dataMatrix = np.full((ncVarData_nrRows, ncVarData_nrCols, NRtimePoints), np.nan)

            timeSeries_dataMatrix[:, :, timePointNR - 1] = ncVarData
            timeSeries_QAMatrix[:, :, timePointNR - 1] = ncQAData

            YRNR = dateYR - min(stats_YRs) + 1
            stats_nrDataPerYR[YRNR] += np.sum(~np.isnan(ncVarData))
            stats_nrDataPerDOY[DOY] += np.sum(~np.isnan(ncVarData))

    # Summing data points spatial
    matrix_nrDataPerPIX = np.sum(~np.isnan(timeSeries_dataMatrix), axis=2)
    matrix_pixWithDataRow, matrix_pixWithDataCol = np.where(matrix_nrDataPerPIX >= 1)
    stats_NRpixWithData = len(matrix_pixWithDataCol)

    # Get full time series data
    timeSeries_dataFull = np.full((stats_NRpixWithData, NRtimePoints), np.nan)
    timeSeries_QAFull = np.full((stats_NRpixWithData, NRtimePoints), np.nan)

    for pixWithDataNR in range(stats_NRpixWithData):
        timeSeries_dataFull[pixWithDataNR, :] = timeSeries_dataMatrix[matrix_pixWithDataRow[pixWithDataNR],
                                                matrix_pixWithDataCol[pixWithDataNR], :]
        timeSeries_QAFull[pixWithDataNR, :] = timeSeries_QAMatrix[matrix_pixWithDataRow[pixWithDataNR],
                                              matrix_pixWithDataCol[pixWithDataNR], :]

    stats_nrDataPerPIX = np.sum(~np.isnan(timeSeries_dataFull), axis=1)

    # Get start and end pixel coordinates for storing
    dataFiles_startPixCoorNRs = np.arange(0, stats_NRpixWithData, 10000)
    dataFiles_endPixCoorNRs = np.unique(np.append(np.arange(10000, stats_NRpixWithData, 10000), stats_NRpixWithData))
    dataFiles_NRdataFiles = len(dataFiles_startPixCoorNRs)

    dataFiles_allFileNames = []

    for dataFileNR in range(dataFiles_NRdataFiles):
        timeSeries_data = timeSeries_dataFull[dataFiles_startPixCoorNRs[dataFileNR]:dataFiles_endPixCoorNRs[dataFileNR],
                          :]
        timeSeries_QA = timeSeries_QAFull[dataFiles_startPixCoorNRs[dataFileNR]:dataFiles_endPixCoorNRs[dataFileNR], :]
        startPixCoorNR = dataFiles_startPixCoorNRs[dataFileNR]
        endPixCoorNR = dataFiles_endPixCoorNRs[dataFileNR]
        fileName_timeSeriesData_pixRange = f"{fileName_timeSeriesData}_pixNRs_{dataFiles_startPixCoorNRs[dataFileNR]}_{dataFiles_endPixCoorNRs[dataFileNR]}"

        sio.savemat(fileName_timeSeriesData_pixRange,
                    {'timeSeries_data': timeSeries_data, 'timeSeries_QA': timeSeries_QA,
                     'startPixCoorNR': startPixCoorNR, 'endPixCoorNR': endPixCoorNR})
        dataFiles_allFileNames.append(fileName_timeSeriesData_pixRange)

    # Save null info
    sio.savemat(fileName_timeSeriesNullInfo,
                {'matrix_nrDataPerPIX': matrix_nrDataPerPIX, 'matrix_pixWithDataRow': matrix_pixWithDataRow,
                 'matrix_pixWithDataCol': matrix_pixWithDataCol, 'timeSeries_dateNR': timeSeries_dateNR,
                 'timeSeries_DOY': timeSeries_DOY, 'stats_nrDataPerPIX': stats_nrDataPerPIX, 'stats_YRs': stats_YRs,
                 'stats_nrDataPerYR': stats_nrDataPerYR, 'stats_nrDataPerDOY': stats_nrDataPerDOY,
                 'stats_NRpixWithData': stats_NRpixWithData, 'dataFiles_allFileNames': dataFiles_allFileNames,
                 'dataFiles_NRdataFiles': dataFiles_NRdataFiles, 'dataFiles_startPixCoorNRs': dataFiles_startPixCoorNRs,
                 'dataFiles_endPixCoorNRs': dataFiles_endPixCoorNRs})

    # Display done
    print(
        f"storeTimeSeriesAndNullInfo_MerisModis done {lakeMetaData['name'][lakeNR]} ({lakeMetaData['short_name'][lakeNR]}): - FULLY")

    return fileName_timeSeriesNullInfo, dataFiles_allFileNames, dataFiles_NRdataFiles