import pandas as pd
import os
from datetime import datetime
from time import time
from functions import storeTimeSeriesAndNullInfo_CCI

def CCISTEP1_storeTimeSeriesAndNullInfoMerisModis_parallelLakes(minLakeListCoor, maxLakeListCoor):
    # Define lake lists
    lakeNRlist_ARRAY = list(range(1, 2025))
    lakeNRlist_MODIS = [1, 2, 5, 11, 14, 19, 22, 25, 36, 39, 40, 41, 50, 52, 56, 57, 64, 65, 68, 72, 74, 91, 110, 111,
                        120, 124, 134, 135, 165, 166, 172, 176, 186, 215, 217, 297, 298, 371, 563, 631, 649, 731, 744,
                        782, 789, 837, 1225, 1230, 1234]
    lakeNRlist_CH = [297, 322, 653, 674, 1335, 1819]
    lakeNRlist_spInterest = [18, 21, 23, 28, 150, 281, 310]

    lakeNRlist_ALL = sorted(set(lakeNRlist_ARRAY + lakeNRlist_MODIS + lakeNRlist_CH + lakeNRlist_spInterest))

    # Load metadata
    file_path = os.path.join('data', 'CCI_lakes_ShpFls', 'lakescci_v2.0.2_data-availability_gQHtOE6.csv')
    lakeMetaData = pd.read_csv(file_path)

    varName = 'chla_mean'  # or 'lwlr_quality_flag'
    skipLakes = False
    skipLakeNRs = []

    minDateNR_FULL = 731335
    maxDateNR_FULL = 738886

    # Loop over lakes
    for lakeListCoor in range(minLakeListCoor, maxLakeListCoor + 1):
        lakeNR = lakeNRlist_ALL[lakeListCoor]

        lake_id = lakeMetaData.loc[lakeNR, 'id']
        lake_short_name = lakeMetaData.loc[lakeNR, 'short_name']
        lake_name = lakeMetaData.loc[lakeNR, 'name']

        print(f"timeSeriesAndNullInfo: {lakeNR} - {lake_id} - {lake_short_name} ({lake_name})")
        print(f"start: {datetime.now().strftime('%d/%m/%y - %H:%M')}")

        if lakeNR not in skipLakeNRs or not skipLakes:
            start_time = time()

            replaceTimeSeriesAndNullInfo = False
            # Call your actual function here
            fileName_timeSeriesNullInfo, dataFiles_allFileNames, dataFiles_NRdataFiles = storeTimeSeriesAndNullInfo_CCI(
                lakeNR, lakeMetaData, varName, minDateNR_FULL, maxDateNR_FULL, replaceTimeSeriesAndNullInfo
            )

            elapsed_time = time() - start_time
            print(f"Processed in {elapsed_time:.2f} seconds")

    print('done')