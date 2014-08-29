#include "eigen.h"

void Eigen::DecrSortEigenStuff(void)
{
    Tridiagonal(); //diagonalize the matrix.
    QLAlgorithm(); //
    DecreasingSort();
    GuaranteeRotation();
}

void Eigen::Tridiagonal(void)
{
    double fM00 = mElement[0][0];
    double fM01 = mElement[0][1];
    double fM02 = mElement[0][2];
    double fM11 = mElement[1][1];
    double fM12 = mElement[1][2];
    double fM22 = mElement[2][2];

    m_afDiag[0] = fM00;
    m_afSubd[2] = 0;
    if (fM02 != (double)0.0)
    {
        double fLength = sqrt(fM01*fM01+fM02*fM02);
        double fInvLength = ((double)1.0)/fLength;
        fM01 *= fInvLength;
        fM02 *= fInvLength;
        double fQ = ((double)2.0)*fM01*fM12+fM02*(fM22-fM11);
        m_afDiag[1] = fM11+fM02*fQ;
        m_afDiag[2] = fM22-fM02*fQ;
        m_afSubd[0] = fLength;
        m_afSubd[1] = fM12-fM01*fQ;
        mElement[0][0] = (double)1.0;
        mElement[0][1] = (double)0.0;
        mElement[0][2] = (double)0.0;
        mElement[1][0] = (double)0.0;
        mElement[1][1] = fM01;
        mElement[1][2] = fM02;
        mElement[2][0] = (double)0.0;
        mElement[2][1] = fM02;
        mElement[2][2] = -fM01;
        m_bIsRotation = false;
    }
    else
    {
        m_afDiag[1] = fM11;
        m_afDiag[2] = fM22;
        m_afSubd[0] = fM01;
        m_afSubd[1] = fM12;
        mElement[0][0] = (double)1.0;
        mElement[0][1] = (double)0.0;
        mElement[0][2] = (double)0.0;
        mElement[1][0] = (double)0.0;
        mElement[1][1] = (double)1.0;
        mElement[1][2] = (double)0.0;
        mElement[2][0] = (double)0.0;
        mElement[2][1] = (double)0.0;
        mElement[2][2] = (double)1.0;
        m_bIsRotation = true;
    }
}

bool Eigen::QLAlgorithm(void)
{
    const int iMaxIter = 32;

    for (int i0 = 0; i0 <3; i0++)
    {
        int i1;
        for (i1 = 0; i1 < iMaxIter; i1++)
        {
            int i2;
            for (i2 = i0; i2 <= (3-2); i2++)
            {
                double fTmp = fabs(m_afDiag[i2]) + fabs(m_afDiag[i2+1]);
                if ( fabs(m_afSubd[i2]) + fTmp == fTmp )
                    break;
            }
            if (i2 == i0)
            {
                break;
            }

            double fG = (m_afDiag[i0+1] - m_afDiag[i0])/(((double)2.0) * m_afSubd[i0]);
            double fR = sqrt(fG*fG+(double)1.0);
            if (fG < (double)0.0)
            {
                fG = m_afDiag[i2]-m_afDiag[i0]+m_afSubd[i0]/(fG-fR);
            }
            else
            {
                fG = m_afDiag[i2]-m_afDiag[i0]+m_afSubd[i0]/(fG+fR);
            }
            double fSin = (double)1.0, fCos = (double)1.0, fP = (double)0.0;
            for (int i3 = i2-1; i3 >= i0; i3--)
            {
                double fF = fSin*m_afSubd[i3];
                double fB = fCos*m_afSubd[i3];
                if (fabs(fF) >= fabs(fG))
                {
                    fCos = fG/fF;
                    fR = sqrt(fCos*fCos+(double)1.0);
                    m_afSubd[i3+1] = fF*fR;
                    fSin = ((double)1.0)/fR;
                    fCos *= fSin;
                }
                else
                {
                    fSin = fF/fG;
                    fR = sqrt(fSin*fSin+(double)1.0);
                    m_afSubd[i3+1] = fG*fR;
                    fCos = ((double)1.0)/fR;
                    fSin *= fCos;
                }
                fG = m_afDiag[i3+1]-fP;
                fR = (m_afDiag[i3]-fG)*fSin+((double)2.0)*fB*fCos;
                fP = fSin*fR;
                m_afDiag[i3+1] = fG+fP;
                fG = fCos*fR-fB;
                for (int i4 = 0; i4 < 3; i4++)
                {
                    fF = mElement[i4][i3+1];
                    mElement[i4][i3+1] = fSin*mElement[i4][i3]+fCos*fF;
                    mElement[i4][i3] = fCos*mElement[i4][i3]-fSin*fF;
                }
            }
            m_afDiag[i0] -= fP;
            m_afSubd[i0] = fG;
            m_afSubd[i2] = (double)0.0;
        }
        if (i1 == iMaxIter)
        {
            return false;
        }
    }
    return true;
}

void Eigen::DecreasingSort(void)
{
    //sort eigenvalues in decreasing order, e[0] >= ... >= e[iSize-1]
    for (int i0 = 0, i1; i0 <= 3-2; i0++)
    {
        // locate maximum eigenvalue
        i1 = i0;
        double fMax = m_afDiag[i1];
        int i2;
        for (i2 = i0+1; i2 < 3; i2++)
        {
            if (m_afDiag[i2] > fMax)
            {
                i1 = i2;
                fMax = m_afDiag[i1];
            }
        }

        if (i1 != i0)
        {
            // swap eigenvalues
            m_afDiag[i1] = m_afDiag[i0];
            m_afDiag[i0] = fMax;
            // swap eigenvectors
            for (i2 = 0; i2 < 3; i2++)
            {
                double fTmp = mElement[i2][i0];
                mElement[i2][i0] = mElement[i2][i1];
                mElement[i2][i1] = fTmp;
                m_bIsRotation = !m_bIsRotation;
            }
        }
    }
}


void Eigen::GuaranteeRotation(void)
{
    if (!m_bIsRotation)
    {
        // change sign on the first column
        for (int iRow = 0; iRow <3; iRow++)
        {
            mElement[iRow][0] = -mElement[iRow][0];
        }
    }
}
