#ifndef EIGEN_H
#define EIGEN_H

#include <cmath>

class Eigen
{
public:
    void DecrSortEigenStuff(void);
    void Tridiagonal(void);
    bool QLAlgorithm(void);
    void DecreasingSort(void);
    void GuaranteeRotation(void);

    double mElement[3][3];
    double m_afDiag[3];
    double m_afSubd[3];
    bool m_bIsRotation;
};

#endif  // EIGEN_H
