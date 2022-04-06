This code corresponds to the paper 'Maximum Likelihood Remission Calibration for Groups of Heterogeneous Laser Scanners' (IROS 2015). If you use this software for scientific publication please cite the paper:

@inproceedings{steder15icra,
  author = {Steder, Bastian and Ruhnke, Michael and K{\"u}mmerle, Rainer and Burgard,
     Wolfram},
  booktitle = {Proc.~of the IEEE Int.~Conf.~on Robotics \& Automation (ICRA)},
  year = {2015},
  title = {Maximum Likelihood Remission Calibration for Groups of Heterogeneous Laser
     Scanners}
} 

-----

How to compile (tested under Ubuntu 14.04, 64Bit):
You will need to have the Point Cloud Library (PCL) installed (http://pointclouds.org/)
- Unpack the code and go into the folder aisRemissionCalibration
- mkdir build
- cd build
- cmake ..
- make

Prepare a dataset:
An example dataset can be downloaded here:
http://www.informatik.uni-freiburg.de/~steder/datasets.html
The format of the point clouds is *.pcd, which is the format of the Point Cloud Library (PCL).
Necessary fields in the point cloud are x, y, z, intensity, laser, which ist the 3D position for each point,
the measured remission value and the laser id. In addition, vp_x, vp_y, vp_z can be used to provide the sensor
position from which the point was observed. If this is not provided, the system uses the VIEWPOINT of the
point cloud as the sensor position for all included points.
The example dataset also includes the field dynProb. This is the probability that a point belongs to a dynamic
object and the system will ignore every point with a value above 0.3.
The only entry necessary in the *_info.dat files is 'SLAM:', which gives the SLAM solution for the dataset
as x, y, z, roll, pitch, yaw (in meters and rad).

Run the code:
- bin/mapCellsExtractor -d <your dataset folder>
  This call extracts datapoints for the optimization from the dataset
- bin/mapCellsExtractor -d <your dataset folder>
- bin/remissionCalibration -d <your dataset folder> -m 4
  This call performs the actual optimization procedure and might take a few minutes.
  The -m parameter defines the mode of the optimization as described in the paper.
  We recommend mode 2. bin/remissionCalibration -h shows a help with more options.
  By default the calibration uses only one cell for the incidence angle, effectively ignoring it.
  This is useful if you only want to calibrate regarding range and laser id, without the need to
  extract normals. Alternatively you can use
- bin/remissionCalibration -d <your dataset folder> -m 4 -y 10
  This performs the optimization including incidence angles.

The result is stored in calibrationValues_mode2.txt.
You can use the following command to write out corrected point clouds:

- bin/remissionCalibration -d <your dataset folder> -n -c calibrationValues_mode2.txt

To visualize a *.pcd you can use 'pcl_viewer', which comes with PCL.

To use calibrated remissions in your code, you have to use the class RemissionCalibrationResult
from src/remissionCalibrationResult.h.

bool readCalibrationValuesFromFile(const std::string& fileName);
reads a calibration file and
float getCalibratedIntensity(float range, float incidenceAngle, float intensity, int laser) const;
calculates the corrected remission for a point.

This software is licenced under the Creative Commons (Attribution-NonCommercial-ShareAlike):
http://creativecommons.org/licenses/by-nc-sa/4.0

Author: Bastian Steder <steder@informatik.uni-freiburg.de>
