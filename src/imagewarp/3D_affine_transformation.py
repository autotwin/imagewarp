#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  : APolonsky
# Created On  : February 9, 2023
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script uses a list of reference points in two datasets to determine the affine transformation
between the two reference frames. This transformation matrix is then used to map the moving data 
into the target reference frame.
using method described by Lenthe et al.: http://dx.doi.org/10.1107/S1600576715009231
reference data contains points in Matrix A
data to be transformed (moving data) contains points in Matrix B
'''
# ---------------------------------------------------------------------------

import numpy as np
import h5py
import pandas as pd

### USER INPUT

### expecting two columns, first with "Reference", second with "Moving" data
matching_features_file = "C:/Users/robomet-team/Desktop/___alignment_test/ListOfMatchingPores.xlsx"
feature_type = "Pore"

## h5 files containing the data
reference_data_file = "C:/Users/robomet-team/Desktop/___alignment_test/ReferenceData/ReferenceData.h5"
moving_data_file = "C:/Users/robomet-team/Desktop/___alignment_test/MovingData/MovingData.h5"

data_file_centroid_ordering = "XYZ" #XYZ or ZYX are acceptable

fill_value = 0

## Array in .h5 file to transform 
## (e.g. "PoreIds" located at "ImageDataContainer/VoxelData/PoreIds" in the .h5 file)
target_voxel_array = "PoreIds"

npy_output_name = "transformed_mapped.npy"

output_vtr = True
if output_vtr:
	from pyevtk.hl import gridToVTK #for writing .vtr files

### END USER INPUT

def get_feature_centroid_matrix(file_path, feature_ids, feature_type):
	with h5py.File(file_path,'r') as file:
		centroids = np.squeeze(file["ImageDataContainer"]["{}Data".format(feature_type)]["Centroids"][:])
		if data_file_centroid_ordering == "ZYX":
			pass
		elif data_file_centroid_ordering == "XYZ":
			centroids[:,[0,2]] = centroids[:,[2,0]] ##swap from XYZ to ZYX ordering
		else:
			raise ValueError("Centroid ordering in h5 file must be XYZ or ZYX")
		feature_centroids = centroids[feature_ids]
		#add column of 1s
		feature_matrix = np.concatenate((feature_centroids, 
										np.expand_dims(np.asarray([1]*feature_centroids.shape[0]), axis=1)),
										axis=1)
		return feature_matrix

def get_transformation_matrix(A,B):
	A_T = np.transpose(A)
	B_T = np.transpose(B)
	T = np.matmul(np.matmul(A_T,B),np.linalg.inv(np.matmul(B_T,B)))
	return T

def get_voxel_resolution(file_path):
	with h5py.File(file_path,'r') as file:
		resolution = np.squeeze(file["ImageDataContainer"]["_Geometry"]["Spacing"][:])
		dx, dy, dz = resolution[0], resolution[1], resolution[2]
		return np.array([dz, dy, dx])

def get_volume_dimensions(file_path):
	with h5py.File(file_path,'r') as file:
		dims = np.squeeze(file["ImageDataContainer"]["_Geometry"]["Dimensions"][:])
		nx, ny, nz = dims[0], dims[1], dims[2]
		return np.array([nz, ny, nx])

def get_image_data(file_path, voxel_array_name):
	with h5py.File(file_path,'r') as file:
		array = np.squeeze(file["ImageDataContainer"]["VoxelData"][voxel_array_name][:])
		return array

def write_vtr(data, res, file_name, data_name):
	### technically swaps X and Z axis but these .vtr files are really just for quick checks
	### Basically the pyevtk package reads in numpy array as XYZ, not ZYX as python does
	### Could swap data array axes first to fix this.
	### Better solution is to put the transformed arrays back into an h5 file in same ZYX ordering
	print("Writing VTR file {}.vtr".format(file_name))
	dx, dy, dz = res[0], res[1], res[2]
	nx, ny, nz = data.shape[0], data.shape[1], data.shape[2]
	z0, y0, x0 = 0,0,0

	# Coordinates
	x = np.arange(x0, (nx+1)*dx, dx, dtype='float64')
	y = np.arange(y0, (ny+1)*dy, dy, dtype='float64')
	z = np.arange(z0, (nz+1)*dz, dz, dtype='float64')

	gridToVTK(file_name, x, y, z, cellData = {data_name: data})

if __name__ == "__main__":
	#get matching point locations from excel file
	df = pd.read_excel(matching_features_file)
	ref_id = df['Reference'].to_numpy()
	mov_id = df['Moving'].to_numpy()

	# format points correctly for the matrix algebra (keeping things ZYX ordering)
	A = get_feature_centroid_matrix(reference_data_file, ref_id, feature_type)
	B = get_feature_centroid_matrix(moving_data_file, mov_id, feature_type)

	#solve for transformation
	T = get_transformation_matrix(A,B)
	T_inv = np.linalg.inv(T)

	# get resolutions of array in dz, dy, dx order as array
	ref_res = get_voxel_resolution(reference_data_file)
	mov_res = get_voxel_resolution(moving_data_file)

	# get resolutions of array in nz, ny, nx order as array
	ref_dim = get_volume_dimensions(reference_data_file)
	mov_dim = get_volume_dimensions(moving_data_file)

	# get image stack to transform
	moving_data = get_image_data(moving_data_file, target_voxel_array)

	### Reference Grid Points Positions (real world coordinates)
	z = np.arange(0, ref_res[0]*ref_dim[0], ref_res[0], dtype=np.int64)
	y = np.arange(0, ref_res[1]*ref_dim[1], ref_res[1], dtype=np.int64)
	x = np.arange(0, ref_res[2]*ref_dim[2], ref_res[2], dtype=np.int64)

	# Fill whole volume with NAN to be filled with mapped data
	transformed_data = np.full((z.shape[0],y.shape[0],x.shape[0]), 0, dtype=int)

	## iterate through every real world position defined by reference grid
	## very slow but uses basically no RAM
	for k in range(0, ref_dim[0]):
		print("Slice {} of {}".format(k+1, ref_dim[0])) #status update
		for j in range(0, ref_dim[1]):
			for i in range(0, ref_dim[2]):
				# find real world position in moving dataset
				mapped_position = np.matmul(T_inv,np.array([z[k],y[j],x[i],1]))[:3]
				# convert to indices in moving dataset (moving dataset resolution)
				mapped_idx = np.rint(np.divide(mapped_position, mov_res)).astype(np.int64)

				# replace non-physical index locations with a user-defined value
				if any(val<0 for val in mapped_idx):
					transformed_data[k,j,i] = fill_value
				elif (mapped_idx[0] >= mov_dim[0]-1) or (mapped_idx[1] >= mov_dim[1]-1) or (mapped_idx[2] >= mov_dim[2]-1):
					transformed_data[k,j,i] = fill_value
				# otherwise copy data in moving dataset at mapped index location
				# to corrseponding index location in the target reference frame
				else:
					transformed_data[k,j,i] = moving_data[mapped_idx[0], mapped_idx[1], mapped_idx[2]]

	np.save(npy_output_name, transformed_data)
	# write out a file to visualize
	if output_vtr:
		write_vtr(transformed_data, ref_res, "transformed_VTR", target_voxel_array)

		ref_data = get_image_data(reference_data_file, target_voxel_array)
		write_vtr(ref_data, ref_res, "reference_VTR", target_voxel_array)