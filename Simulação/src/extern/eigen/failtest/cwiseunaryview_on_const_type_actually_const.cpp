#include "../Eigen/Core"

#ifdef EIGEN_SHOULD_FAIL_TO_BUILD
#define CV_QUALIFIER const
#else
#define CV_QUALIFIER
#endif

using namespace Eigen;

void foo() {
  MatrixXf m;
  CwiseUnaryView<internal::scalar_real_ref_op<float>, CV_QUALIFIER MatrixXf>(m).coeffRef(0, 0) = 1.0f;
}

int main() {}
