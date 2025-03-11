# AAE6102 Assignment 1: Satellite Communication and Navigation

**Semester:** 2024/25 Semester 2

## Task 1 – Acquisition

**Objective:** Process IF data using a GNSS SDR and generate initial acquisition results.  The result file is in `Task1`.

**Table 1. Initial Acquisition Results of Open-Sky Dataset**

| SV | SNR       | Doppler | CodeDelay | FineFreq    |
|----|-----------|---------|-----------|-------------|
| 3  | 18.0953   | 1000    | 3683      | 4580990     |
| 4  | 17.2936   | -3000   | 12701     | 4576905     |
| 16 | 26.4349   | 0       | 26051     | 4579695     |
| 22 | 19.8321   | 1500    | 2610      | 4581565     |
| 26 | 27.2060   | 2000    | 57908     | 4581835     |
| 27 | 22.7206   | -3000   | 49778     | 4576775     |
| 31 | 24.4017   | 1000    | 39064     | 4581045     |
| 32 | 22.1983   | 3500    | 20170     | 4583345     |

**Description:** Open-Sky: SV (Satellite IDs): [3, 4, 16, 22, 26, 27, 31, 32]. SNR (Signal-to-Noise Ratio): [18.0953, 17.2936, 26.4349, 19.8321, 27.2060, 22.7206, 24.4017, 22.1983]. Doppler: [1000, -3000, 0, 1500, 2000, -3000, 1000, 3500]. CodeDelay: [3683, 12701, 26051, 2610, 57908, 49778, 39064, 20170]. FineFreq: [4580990, 4576905, 4579695, 4581565, 4581835, 4576775, 4581045, 4583345].


**Table 2. Initial Acquisition Results of Urban Dataset**

| SV | SNR       | Doppler | CodeDelay | FineFreq |
|----|-----------|---------|-----------|----------|
| 1  | 42.6344   | 1000    | 22742     | 1200     |
| 3  | 29.3771   | 4500    | 1154      | 4285     |
| 7  | 19.9028   | 500     | 10811     | 365      |
| 11 | 23.0521   | 500     | 24851     | 405      |
| 22 | 17.7927   | 3500    | 2050      | 3315     |

**Description:** Urban: SV (Satellite IDs): [1, 3, 7, 11, 22]. SNR (Signal-to-Noise Ratio): [42.6344, 29.3771, 19.9028, 23.0521, 17.7927]. Doppler: [1000, 4500, 500, 500, 3500]. CodeDelay: [22742, 1154, 10811, 24851, 2050]. FineFreq: [1200, 4285, 365, 405, 3315].


## Task 2 – Tracking

**Objective:** Adapt the tracking loop (DLL) to generate correlation plots and analyze tracking performance. Discuss the impact of urban interference on correlation peaks. (Multiple correlators must be implemented.) The result file is in `Task2`.

**(Include Figures 1 and 2 from the PDF here.  You'll need to save these images separately and reference them in your markdown using the standard Markdown image syntax.)**

**Discussion on Urban Interference:**

### Multipath Effects:
In urban environments, signals can bounce off buildings and other structures, causing multipath interference. This results in multiple signal paths arriving at different times, which can lead to smearing of correlation peaks.

### Tracking Difficulties:
Under significant urban interference, the DLL may have difficulty locking onto the signal due to competing signals, leading to a lower correlation peak and wider spread.

### Mitigation Strategies:
* **Adaptive Algorithms:** Consider using adaptive filtering techniques to improve performance in multipath environments.
* **Signal Processing Enhancements:** Incorporate additional techniques like signal averaging or combining measurements from multiple satellites to improve robustness against interference.


## Task 3 – Navigation Data Decoding

**Objective:** Decode the navigation message and extract key parameters, such as ephemeris data, for at least one satellite. The result file is in `Task3`.

* **`eph` (Ephemeris):** Represents the satellite's ephemeris data. Contains orbital information, crucial for accurate positioning.
* **`sbf` (Subframe):** Represents subframe data from the navigation message. Contains information about the satellite's status, health, timing, and other navigation-related details.
* **`TckResult_Eph` (Tracking Result Ephemeris):** Stores tracking results (signal strength, frequency offset, code phase, etc.), used for further navigation data processing.  *(Note: This file is too large to upload.)*


## Task 4 – Position and Velocity Estimation

**Objective:** Using pseudorange measurements from tracking, implement the Weighted Least Squares (WLS) algorithm to compute the user's position and velocity. Plot the user position and velocity. Compare the results with the ground truth. Discuss the impact of multipath effects on the WLS solution.

**(Include Figures 3 and 4 from the PDF here. You'll need to save these images separately and reference them in your Markdown using the standard Markdown image syntax.)**

**Effects of Multipath on Weighted Least Squares (WLS) Solutions**

1. **Increased Measurement Errors:** Multipath causes time delays in received signals, leading to pseudorange errors and inaccurate WLS estimates.
2. **Impact on Weight Selection:** Multipath increases variance in some measurements, affecting weight assignments and potentially reducing accuracy.
3. **Poor Model Fitting:** Multipath-induced data deviations lead to poorer model fits and less accurate solutions.
4. **Impact on Estimated Results:** Multipath causes coordinate estimates to deviate from true values, especially in obstructed environments.

**Some Methods to Mitigate Multipath Effects**

1. **Signal Processing Techniques:** Use advanced techniques like multi-channel receivers or differential GNSS.
2. **Use Appropriate Weights:** Dynamically adjust weights based on measurement quality.
3. **Improve Models:** Develop models that incorporate multipath errors.
4. **Environmental Adaptation:** Choose suitable locations and use multipath-resistant antennas.


## Task 5 – Kalman Filter-Based Positioning

**Objective:** Develop an Extended Kalman Filter (EKF) using pseudorange and Doppler measurements to estimate user position and velocity.

