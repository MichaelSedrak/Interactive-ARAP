#pragma once
#include "ProcrustesAligner.h"

Matrix4f ProcrustesAligner::estimatePose(const std::vector<Vector3f>& sourcePoints, const std::vector<Vector3f>& targetPoints) {
	ASSERT(sourcePoints.size() == targetPoints.size() && "The number of source and target points should be the same, since every source point is matched with corresponding target point.");

	// We estimate the pose between source and target points using Procrustes algorithm.
	// Our shapes have the same scale, therefore we don't estimate scale. We estimated rotation and translation
	// from source points to target points.

	auto sourceMean = computeMean(sourcePoints);
	auto targetMean = computeMean(targetPoints);

	std::cout << "Source/ Target mean computed." << std::endl;

	Matrix3f rotation = estimateRotation(sourcePoints, sourceMean, targetPoints, targetMean);
	Vector3f translation = computeTranslation(sourceMean, targetMean, rotation);

	// TODO: Compute the transformation matrix by using the computed rotation and translation.
	// You can access parts of the matrix with .block(start_row, start_col, num_rows, num_cols) = elements
	Matrix4f estimatedPose;
	estimatedPose.block<3, 3>(0, 0) = rotation;
	estimatedPose.block<1, 4>(3, 0) << 0, 0, 0, 1;
	estimatedPose.block<3, 1>(0, 3) = translation;
	return estimatedPose;
}

Vector3f ProcrustesAligner::computeMean(const std::vector<Vector3f>& points) {
	// TODO: Compute the mean of input points.

	Vector3f temp(0.0, 0.0, 0.0);
	int i = 0;

	for (auto& point : points) {
		temp += point;
		i++;
	}

	return temp /= i;
}

Matrix3f ProcrustesAligner::estimateRotation(const std::vector<Vector3f>& sourcePoints, const Vector3f& sourceMean, const std::vector<Vector3f>& targetPoints, const Vector3f& targetMean) {
	// TODO: Estimate the rotation from source to target points, following the Procrustes algorithm. 
	// To compute the singular value decomposition you can use JacobiSVD() from Eigen.

	// Zero center data and concatenate to matrix
	// Target: x, source: x_hat

	MatrixXf source(4, 3);
	MatrixXf target(4, 3);

	int i = 0;
	int j = 0;

	for (auto& point : sourcePoints) {
		source.row(i) = (point - sourceMean);
		i++;
	}

	for (auto& point : targetPoints) {
		target.row(j) = (point - targetMean);
		j++;
	}

	std::cout << "Source - mean" << std::endl;
	std::cout << sourceMean << std::endl;
	std::cout << source << std::endl;
	std::cout << "Target - mean" << std::endl;
	std::cout << targetMean << std::endl;
	std::cout << target << std::endl;

	// Calculate Cross-Covariance Matrix
	//x^{T}*x_hat

	JacobiSVD<MatrixXf> svd(target.transpose() * source, ComputeFullU | ComputeFullV);
	std::cout << "U" << std::endl;
	std::cout << svd.matrixU() << std::endl;
	std::cout << "V" << std::endl;
	std::cout << svd.matrixV() << std::endl;
	std::cout << "R" << std::endl;
	std::cout << svd.matrixU() * svd.matrixV().transpose() << std::endl;
	std::cout << "R.determinant()" << std::endl;
	std::cout << (svd.matrixU() * svd.matrixV().transpose()).determinant() << std::endl;
	return (svd.matrixU() * svd.matrixV().transpose());
}

Vector3f ProcrustesAligner::computeTranslation(const Vector3f& sourceMean, const Vector3f& targetMean, const Matrix3f& rotation) {
	// TODO: Compute the translation vector from source to target opints.

	return  (-rotation * sourceMean) + targetMean;
}