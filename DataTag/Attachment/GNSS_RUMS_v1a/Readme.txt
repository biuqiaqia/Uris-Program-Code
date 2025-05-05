GNSS-RUMS Code Introduction

Platform: MATLAB

Code environment setup: Unzip "MATLAB_library.rar" and set path to those folders

Main code: "GNSS_RUMS_main_ver1a_demo.m"

Simulation Input: 
1. Agents' true trajectories (Line 20)
	- Rows for different agents
		- Column 1: 	Index of epoch
		- Column 2-4:	Agent true position in latitude, longitude, ellipsoid height [degree,degree,meter]
		- Column 5-7:   Agent true velocity in East-North-Up (ENU) coordinate [m/s,m/s,m/s]
2. Ephemeris data direction (Line 21): Point to the ephermeris data in Rinex 3.xx format of GPS, other constellations in the same folder will be loaded automatically.
3. GPS time, GPS week (Line 22): Week number
4. GPS time, GPS seconds of week (Line 23): Different times of week [seconds], aligned with the agent true trajectory data
5. 3D building model direction (Line 24): Point to the model file in Google Earth KML format
* currently only support GPS satellites and Beidou C01-C37 satellites.
* default is runing with MATLAB parallel pool, can be changed to single processor by changing Line 81 "parfor" -> "for".
* carrier-phase simulation is not supported and still under developing.

Simulation Output "Sim_data": 
	- Column 1: GPS time in [GPS week, Seconds of week]
	- Column 2 to end: Simulated measurements for different agents
	- Rows: Simulated measurements for different epochs
	- Inside:	- Column 1: Satellite index, 1-32 for GPS G01-G32 and 87-123 for Beidou C01-C37
			- Column 2: Satellite pseudorange measurements (meter) with simulated errors
			- Column 3: Satellite measurement carrier-to-noise ratio (dB-Hz)
			- Column 4: Satellite Doppler shift measurements (Hz)

Intermediate variables during simulation can be found in Line 480-490

Simulation algorithms are based on:
Zhang, G., Xu, B., Ng, H.-F., & Hsu, L.-T. (2021). GNSS RUMS: GNSS Realistic Urban Multiagent Simulator for Collaborative Positioning Research. Remote Sensing, 13(4), 544. 