#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  : APolonsky
# Created On  : March 8, 2023
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script uses a list of feature control points in two datasets to determine the affine transformation
between the two reference frames. This transformation matrix is then used to map the moving data 
into the target reference frame. Works on two sets of image stacks.
using method described by Lenthe et al.: http://dx.doi.org/10.1107/S1600576715009231
reference data contains points in Matrix A
data to be transformed (moving data) contains points in Matrix B
Feature control points are expected in an excel file with column names of 
"PosX", "PosY", and "PosZ" to denote the real-world locations (not index locations)
of control points. Identical features are paired row-wise in both control point files
'''
# ---------------------------------------------------------------------------

import numpy as np
import pandas as pd
import glob
from PIL import Image
import math

### USER INPUT

# reference data (frame we want to map into) A in the paper
ref_points_file = "D:/Al2219_DCT-TriBeam/BSE-CT_Alignment/CT_Points.xlsx" #locations in real space
ref_image_dir = "D:/Al2219_DCT-TriBeam/BSE-CT_Alignment/CT_Grayscale/"
ref_dx, ref_dy, ref_dz = 1.5, 1.5, 1.5
ref_file_type = ".tif"
# index bounds in reference dataframe to speed up mapping (use 0,-1 for full range)
start_z, end_z = 150, 500
start_y, end_y = 240, 240+250
start_x, end_x = 315, 315+250

# moving data (frame we want to transform) B in the paper
mov_points_file = "D:/Al2219_DCT-TriBeam/BSE-CT_Alignment/BSE_Points.xlsx" #locations in real space
# mov_image_dir = "D:/Al2219_DCT-TriBeam/BSE-CT_Alignment/BSE_Grayscale/"
# mov_dx, mov_dy, mov_dz = 1600/6144, 1600/6144, 2.0
# mov_file_type = ".tif"
mov_image_dir = "D:/Al2219_DCT-TriBeam/BSE-CT_Alignment/EBSD_Vector/"
mov_dx, mov_dy, mov_dz = 2.0, 2.0, 2.0
mov_file_type = ".tif"

fill_value = 0

save_dir = "D:/Al2219_DCT-TriBeam/BSE-CT_Alignment/"
npy_output_name = "Vector_mapped.npy"

### END USER INPUT

def transform_list_of_points(transformation_matrix, ZYX_points, direction):
	## real space coordinate transformation
	T = transformation_matrix
	T_inv = np.linalg.inv(T)
	transformed_points = []
	for point in ZYX_points:
		if direction == "moving_to_reference":
			### using T and not inverse T because going from moving into reference, not reference into moving
			mapped_point = np.matmul(T,np.array([point[0], point[1], point[2],1]))
		elif direction == "reference_to_moving":
			mapped_point = np.matmul(T_inv,np.array([point[0], point[1], point[2],1]))
		transformed_points.append(mapped_point)
	transformed_points = np.asarray(transformed_points)
	return transformed_points[:,:3]


def get_feature_centroid_matrix_XYZ(file_path):

	df = pd.read_excel(file_path)
	x_pos = df['PosX'].to_numpy()
	y_pos = df['PosY'].to_numpy()
	z_pos = df['PosZ'].to_numpy()

	centroids = np.stack((x_pos,y_pos,z_pos), axis=-1) #X,Y,Z ordering for math

	#add column of 1s
	feature_matrix = np.concatenate((centroids, 
									np.expand_dims(np.asarray([1]*centroids.shape[0]), axis=1)),
									axis=1)
	return feature_matrix

def get_transformation_matrix(A,B):
	A_T = np.transpose(A)
	B_T = np.transpose(B)
	T = np.matmul(np.matmul(A_T,B),np.linalg.inv(np.matmul(B_T,B)))
	return T

if __name__ == "__main__":
	# format points correctly for the matrix algebra (keeping things ZYX ordering)
	A = get_feature_centroid_matrix_XYZ(ref_points_file)
	B = get_feature_centroid_matrix_XYZ(mov_points_file)

	#solve for transformation
	T = get_transformation_matrix(A,B)
	T_inv = np.linalg.inv(T)

	## Decompositon of the affine transformation matrix
	# https://colab.research.google.com/drive/1ImBB-N6P9zlNMCBH9evHD6tjk0dzvy1_#scrollTo=sA10qtzcCppr
	#https://math.stackexchange.com/questions/237369/given-this-transformation-matrix-how-do-i-decompose-it-into-translation-rotati
	H = np.copy(T)
	print(H)

	### TRANSLATION
	# Translation (T) from last column
	T = np.eye(4)
	T[:3,3] = H[:3,3]
	# print(T)
	# and linear component (L) from 3x3 upper-left block
	L = H.copy()
	L[:3,3] = 0
	# print(L)
	# check that H = TL
	if not np.allclose(H, T @ L):
		raise ValueError("Error in rotation matrix decomposition")

	### ROTATION
	# calc polar decomposition of L to obtain rotation (R) and stretch (K) matrices
	# L = RK
	# H = TRK
	from scipy.linalg import polar
	R, K = polar(L)
	# The determinant of a rotation matrix must be positive.
	# When the determinant of  𝑅  is negative it is necessary to 
	# change the sign of the linear part of both matrices:
	if np.linalg.det(R) < 0:
		print("Changing sign of L")
		R[:3,:3] = -R[:3,:3]
		K[:3,:3] = -K[:3,:3]
	# print(R)
	# print(K)
	#Check that  𝐿=𝑅𝐾  and  𝐻=𝑇𝑅𝐾 :
	if not np.allclose(L, R @ K):
		raise ValueError("Error in rotation matrix decomposition")
	if not np.allclose(H, T @ R @ K):
		raise ValueError("Error in rotation matrix decomposition")
	### ROTATION CONVERTER
	# https://www.andre-gaschler.com/rotationconverter/
	# https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/
	angle = math.acos((R[0][0]+R[1][1]+R[2][2]-1)/2) * 180/math.pi
	denominator = math.sqrt((R[2][1]-R[1][2])**2 + (R[0][2]-R[2][0])**2 + (R[1][0]-R[0][1])**2)
	i = (R[2][1]-R[1][2])/denominator
	j = (R[0][2]-R[2][0])/denominator
	k = (R[1][0]-R[0][1])/denominator
	print("Angle: {} degrees".format(angle))
	print("Axis [i,j,k] = [{}, {}, {}]".format(i,j,k))
	print(R)

	### SCALE
	# The stretch matrix  𝐾  can be further analyzed to obtain up to three scale matrices.
	# 𝐾=𝑆1𝑆2𝑆3 
	# 𝐻=𝑇𝑅𝑆1𝑆2𝑆3 
	# Each eigenvalue of  𝐾  represents a scale factor  𝑓𝑖 , and its corresponding eigenvector  𝑋𝑖  represents the scaling axis. 
	# Each pair ( 𝑓𝑖,𝑋𝑖 ) can be used to construct a scale matrix:
	# 𝑆𝑖=𝐼+𝑋𝑖𝑋𝑇𝑖(𝑓𝑖−1) 
	# When  𝑓𝑖=1  the result is an identity matrix, so the pair can be discarded.
	f, X = np.linalg.eig(K)
	# print(f)
	# print(X)
	S = []
	for factor, axis in zip(f, X.T):
		if not np.isclose(factor, 1):
		    scale = np.eye(4) + np.outer(axis, axis) * (factor-1)
		    S.append(scale)
		    # print(scale)
    # Check that  𝐾=𝑆1𝑆2𝑆3  and  𝐻=𝑇𝑅𝑆1𝑆2𝑆3 :
	if not np.allclose(K, S[0] @ S[1] @ S[2]):
		raise ValueError("Error in rotation matrix decomposition")
	if not np.allclose(H, T @ R @ S[0] @ S[1] @ S[2]):
		raise ValueError("Error in rotation matrix decomposition")
