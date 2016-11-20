#include "d_vanderpol_linear.h"
#include "src/math/constants.h"

namespace Tasks
{

namespace Discrete
{

using Math::Rand::gaussianVector;


VanDerPolLinear::VanDerPolLinear()
    : DiscreteTask()
    , m_omega(0.1 * Math::Const::PI)
    , m_alpha(1.0)
    , m_beta(0.25)
{
    m_info->setName("Осциллятор Ван-дер-Поля");
    m_info->setType("Л-");

    m_dimY = 2;

    m_dimX      = 2;
    m_meanX0    = Vector(m_dimX);
    m_meanX0[0] = 10.0;
    m_meanX0[1] = -3.0;

    m_dimV  = 2;
    m_meanV = Vector::Zero(m_dimV);

    m_dimW     = 2;
    m_meanW    = Vector(m_dimW);
    m_meanW[0] = 1.0;
    m_meanW[1] = 1.5;

    m_varX0 = Matrix::Zero(m_dimX, m_dimX);
    m_varX0(0, 0) = 5.0;
    m_varX0(1, 1) = 5.0;

    m_varV = Matrix::Identity(m_dimV, m_dimV);

    m_varW = Matrix::Identity(m_dimW, m_dimW);
    m_varW(0, 0) = 4.0;
    m_varW(1, 1) = 4.0;

    (*m_params)["Omega"] = m_omega;
    (*m_params)["Alpha"] = m_alpha;
    (*m_params)["Beta"]  = m_beta;
}

void VanDerPolLinear::loadParams()
{
    m_omega = m_params->at("Omega");
    m_alpha = m_params->at("Alpha");
    m_beta  = m_params->at("Beta");
}

Vector VanDerPolLinear::a(const Vector &x) const
{
    Vector res(m_dimX);

    res[0] = x[1];
    res[1] = -m_omega * m_omega * x[0] + m_alpha * x[1] - m_alpha * m_beta * x[0] * x[0] * x[1];

    return res;
}

Matrix VanDerPolLinear::B(const Vector &x) const
{
    Matrix res = Matrix::Zero(m_dimX, m_dimV);

    res(1, 1) = x[0];

    return res;
}

Vector VanDerPolLinear::tau(const Vector &m, const Matrix & /*D*/) const
{
    Vector res(m_dimX);

    res[0] = m[1];
    res[1] = -m_omega * m_omega * m[0] + m_alpha * m[1] - m_alpha * m_beta * m[0] * m[0] * m[1];

    return res;
}

Matrix VanDerPolLinear::Theta(const Vector &m, const Matrix & /*D*/) const
{
    Matrix b = B(m);

    return b * b.transpose();
}

Matrix VanDerPolLinear::A(const Vector &m, const Matrix & /*D*/) const
{
    Matrix res(m_dimX, m_dimX);

    res(0, 0) = 0.0;
    res(0, 1) = 1.0;
    res(1, 0) = -m_omega * m_omega - 2.0 * m_alpha * m_beta * m[0] * m[1];
    res(1, 1) = m_alpha * (1 - m_beta * m[0] * m[0]);

    return res;
}


Vector VanDerPolLinear::c(const Vector &x) const
{
    return x + gaussianVector(m_meanW, m_varW);
}

Vector VanDerPolLinear::h(const Vector &m, const Matrix & /* D*/) const
{
    return m + m_meanW;
}

Matrix VanDerPolLinear::G(const Vector & /*m*/, const Matrix & /*D*/) const
{
    Matrix res = Matrix::Zero(2, 2); // WARNING (size = ?)

    res(0, 0) = 1.0;

    return res;
}

Matrix VanDerPolLinear::F(const Vector & /*m*/, const Matrix &D) const
{
    return D + m_varW;
}

Matrix VanDerPolLinear::dadx(const Vector &x) const
{
//    double e = exp(-BB * x[2]);
//    Matrix res(m_dimX, m_dimX);

//    res(0, 0) = 1.0 - 2.0 * m_step * CC * x[0] * e;
//    res(0, 1) = -m_step * GG * cos(x[1]);
//    res(0, 2) = m_step * BB * CC * x[0] * x[0] * e;

//    res(1, 0) = m_step * (CC * k(m_time) * e + (GG / (x[0] * x[0]) + 1.0 / (RR + x[2])) * cos(x[1]));
//    res(1, 1) = 1.0 + m_step * ((GG / x[0] - x[0] / (RR + x[2])) * sin(x[1]));
//    res(1, 2) = m_step * (-CC * BB * k(m_time) * x[0] * e - x[0] * cos(x[1]) / ((RR + x[2]) * (RR + x[2])));

//    res(2, 0) = m_step * sin(x[1]);
//    res(2, 1) = m_step * x[0] * cos(x[1]);
//    res(2, 2) = 1.0;

//    return res;
}

Matrix VanDerPolLinear::dadv(const Vector & /*x*/) const
{
//    return Matrix::Zero(m_dimX, m_dimV);
}

Matrix VanDerPolLinear::dbdx(const Vector &x) const
{
//    double e0 = CC * x[0] * x[0] * exp(-BB * x[2]) * (cos(x[1]) - k(m_time) * sin(x[1]));
//    double e1 = CC * x[0] * x[0] * exp(-BB * x[2]) * (sin(x[1]) - k(m_time) * cos(x[1]));
//    Matrix res(m_dimY, m_dimX);

//    res(0, 0) = 2.0 * e0 * (1.0 + m_meanW[0]) / x[0];
//    res(0, 1) = -e1 * (1.0 + m_meanW[0]);
//    res(0, 2) = -BB * e0 * (1.0 + m_meanW[0]);

//    res(1, 0) = 2.0 * e1 * (1.0 + m_meanW[1]) / x[0];
//    res(1, 1) = e0 * (1.0 + m_meanW[1]);
//    res(1, 2) = -BB * e1 * (1.0 + m_meanW[1]);

//    return res;
}

Matrix VanDerPolLinear::dbdw(const Vector &x) const
{
//    Matrix res(m_dimY, m_dimW);

//    res(0, 0) = CC * x[0] * x[0] * exp(-BB * x[2]) * (cos(x[1]) - k(m_time) * sin(x[1]));
//    res(0, 1) = 0.0;
//    res(0, 2) = 1.0;
//    res(0, 3) = 0.0;

//    res(1, 0) = 0.0;
//    res(1, 1) = CC * x[0] * x[0] * exp(-BB * x[2]) * (sin(x[1]) - k(m_time) * cos(x[1]));
//    res(1, 2) = 0.0;
//    res(1, 3) = 1.0;

//    return res;
}

Vector VanDerPolLinear::b(const Vector &x) const
{
//    double e = exp(-BB * x[2]);
//    Vector w = gaussianVector(m_meanW, m_varW);
//    Vector res(m_dimY);

//    res[0] = CC * (w[0] + 1.0) * x[0] * x[0] * e * (cos(x[1]) - k(m_time) * sin(x[1])) + w[2];
//    res[1] = CC * (w[1] + 1.0) * x[0] * x[0] * e * (sin(x[1]) - k(m_time) * cos(x[1])) + w[3];

//    return res;
}

} // end Tasks::Discrete

} // end Tasks
