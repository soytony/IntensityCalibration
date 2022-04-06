template <typename real, int dimension>
VectorAverage<real, dimension>::VectorAverage()
{
  reset();
}

template <typename real, int dimension>
inline void VectorAverage<real, dimension>::reset()
{
  noOfSamples_ = 0;
  accumulatedWeight_ = 0.0;
  mean_.fill(0);
  covariance_.fill(0);
}

template <typename real, int dimension>
inline void VectorAverage<real, dimension>::add(const Eigen::Matrix<real, dimension, 1>& sample, real weight) {
  if (weight == 0.0f)
    return;
  
  ++noOfSamples_;
  accumulatedWeight_ += weight;
  real alpha = weight/accumulatedWeight_;
  
  Eigen::Matrix<real, dimension, 1> diff = sample - mean_;
  covariance_ = (1.0-alpha)*(covariance_ + alpha * (diff * diff.transpose()));
  
  mean_ += alpha*(diff);

  //if (pcl_isnan(covariance_(0,0)))
  //{
    //cout << PVARN(weight);
    //exit(0);
  //}
}

template <typename real, int dimension>
inline void VectorAverage<real, dimension>::doPCA(Eigen::Matrix<real, dimension, 1>& eigen_values,
                                Eigen::Matrix<real, dimension, 1>& eigen_vector1, Eigen::Matrix<real, dimension, 1>& eigen_vector2) const
{
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, dimension, dimension> > ei_symm(covariance_);
  eigen_values = ei_symm.eigenvalues();
  Eigen::Matrix<real, dimension, dimension> eigen_vectors = ei_symm.eigenvectors();
  
  eigen_vector1 = eigen_vectors.col(0);
  eigen_vector2 = eigen_vectors.col(1);
}

template <typename real, int dimension>
inline void VectorAverage<real, dimension>::doPCA(Eigen::Matrix<real, dimension, 1>& eigen_values, Eigen::Matrix<real, dimension, 1>& eigen_vector1,
                                                  Eigen::Matrix<real, dimension, 1>& eigen_vector2, Eigen::Matrix<real, dimension, 1>& eigen_vector3) const
{
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, dimension, dimension> > ei_symm(covariance_);
  eigen_values = ei_symm.eigenvalues();
  Eigen::Matrix<real, dimension, dimension> eigen_vectors = ei_symm.eigenvectors();
  
  eigen_vector1 = eigen_vectors.col(0);
  eigen_vector2 = eigen_vectors.col(1);
  eigen_vector3 = eigen_vectors.col(2);
}

template <typename real, int dimension>
inline void VectorAverage<real, dimension>::doPCA(Eigen::Matrix<real, dimension, 1>& eigen_values) const
{
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, dimension, dimension> > ei_symm(covariance_, false);
  //Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, dimension, dimension> > ei_symm(covariance_);
  eigen_values = ei_symm.eigenvalues();
}

template <typename real, int dimension>
inline void VectorAverage<real, dimension>::getEigenVector1(Eigen::Matrix<real, dimension, 1>& eigen_vector1) const
{
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, dimension, dimension> > ei_symm(covariance_);
  Eigen::Matrix<real, dimension, dimension> eigen_vectors = ei_symm.eigenvectors();
  eigen_vector1 = eigen_vectors.col(0);
}

#if 0  // Fuck... Does not work always...
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Special cases for real=float & dimension=3 -> Partial specialization does not work with class templates. :( //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
// float //
///////////
template <>
inline void VectorAverage<float, 3>::doPCA(Eigen::Matrix<float, 3, 1>& eigen_values, Eigen::Matrix<float, 3, 1>& eigen_vector1,
                                          Eigen::Matrix<float, 3, 1>& eigen_vector2, Eigen::Matrix<float, 3, 1>& eigen_vector3) const
{
  //cout << "Using specialized 3x3 version of doPCA!\n";
  Eigen::Matrix<float, 3, 3> eigen_vectors;
  eigen33(covariance_, eigen_vectors, eigen_values);
  eigen_vector1 = eigen_vectors.col(0);
  eigen_vector2 = eigen_vectors.col(1);
  eigen_vector3 = eigen_vectors.col(2);
}
template <>
inline void VectorAverage<float, 3>::doPCA(Eigen::Matrix<float, 3, 1>& eigen_values) const
{
  //cout << "Using specialized 3x3 version of doPCA!\n";
  computeRoots (covariance_, eigen_values);
}
template <>
inline void VectorAverage<float, 3>::getEigenVector1(Eigen::Matrix<float, 3, 1>& eigen_vector1) const
{
  //cout << "Using specialized 3x3 version of doPCA!\n";
  Eigen::Matrix<float, 3, 1> eigen_values;
  Eigen::Matrix<float, 3, 3> eigen_vectors;
  eigen33(covariance_, eigen_vectors, eigen_values);
  eigen_vector1 = eigen_vectors.col(0);
}

////////////
// double //
////////////
template <>
inline void VectorAverage<double, 3>::doPCA(Eigen::Matrix<double, 3, 1>& eigen_values, Eigen::Matrix<double, 3, 1>& eigen_vector1,
                                          Eigen::Matrix<double, 3, 1>& eigen_vector2, Eigen::Matrix<double, 3, 1>& eigen_vector3) const
{
  //cout << "Using specialized 3x3 version of doPCA!\n";
  Eigen::Matrix<double, 3, 3> eigen_vectors;
  eigen33(covariance_, eigen_vectors, eigen_values);
  eigen_vector1 = eigen_vectors.col(0);
  eigen_vector2 = eigen_vectors.col(1);
  eigen_vector3 = eigen_vectors.col(2);
}
template <>
inline void VectorAverage<double, 3>::doPCA(Eigen::Matrix<double, 3, 1>& eigen_values) const
{
  //cout << "Using specialized 3x3 version of doPCA!\n";
  computeRoots (covariance_, eigen_values);
}
template <>
inline void VectorAverage<double, 3>::getEigenVector1(Eigen::Matrix<double, 3, 1>& eigen_vector1) const
{
  //cout << "Using specialized 3x3 version of doPCA!\n";
  Eigen::Matrix<double, 3, 1> eigen_values;
  Eigen::Matrix<double, 3, 3> eigen_vectors;
  eigen33(covariance_, eigen_vectors, eigen_values);
  eigen_vector1 = eigen_vectors.col(0);
}
#endif
