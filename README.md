# Uris-Program-Code
Storing the source codes Uris porgram

## DataProcess
### Code RINEX_processing_main_f.m 
is a Matlab code for processing GNSS data and outputting a csv file which contains 
% GPS_Time
% PRN
% Elevation
% Azimuth
% C/N0
% Pr_Residuals
% WLS_Latitude
% WLS_Longitude
% WLS_Height
% GT_Latitude
% GT_Longitude
% GT_Height
% doppler
% SaVX
% SaVY
% SaVZ
% lambda
% corrected_Pr
% SaX
% SaY
% SaZ
% Max_PRN
% Max_Residual
### Folders "goGPS", "fsaxen-ParforProgMon-f99c6f0", "gtsam_toolbox-master"
are all the required environmental configuration files, such as the required functions or toolkits.

## DataTag
### Code zzy_skymask_generation_f.m
is a Matlab code for marking the GNSS data type. Using kml files containing corresponding building information.
type 0 is NLOS
type 1 is LOS
type 2 is NLOS_Refl (NLOS with reflection)
type 3 is NLOS_Diff (NLOS with diffraction)
type 4 is LOS (LOS with diffraction)
Besides, it will also output the distance and point of the reflection or diffraction and diffraction coefficient.
There is also a code that determines whether the signal is obscured through azimuth and elevation angles. It will output "Visible" or "Blocked"
### Folder "Attachment" 
is the required environmental configuration files, such as the required functions or toolkits.

## MLP
### Code MLP_1.py
is a python code. A MLP model for signal classification detached from building data.
This is a preliminary model. Although it has undergone multiple adjustments and improvements, it has not yet been officially confirmed for use.

## Visualized_Path
### Code to kml.f.py
is a python code. Generate the Google earth visual path through the positions of the sampling points and the satellites. Output kml files.

## TagDatabase
### File TagDatabase.csv
is all the labeled GNSS data accumulated as of May 1st.
