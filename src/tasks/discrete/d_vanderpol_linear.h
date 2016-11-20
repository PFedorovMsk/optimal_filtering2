#ifndef DISCRETE_VANDERPOL_LINEAR_H
#define DISCRETE_VANDERPOL_LINEAR_H

#include "src/core/discrete_task.h"
#include "src/math/math.h"

namespace Tasks
{

namespace Discrete
{


/*!
 * \brief Осциллятор Ван-дер-Поля (дискретный, линеаризованный) для фильтров оптимальной структуры.
 */


class VanDerPolLinear : public Core::DiscreteTask
{

public:
    //! \brief Конструктор.
    VanDerPolLinear();

    Vector a(const Vector &x) const override;
    Matrix B(const Vector &x) const override;
    Vector c(const Vector &x) const override;
    Vector b(const Vector &x) const override;
    Vector tau(const Vector &m, const Matrix &D) const override;
    Matrix Theta(const Vector &m, const Matrix &D) const override;
    Matrix A(const Vector &m, const Matrix &D) const override;
    Vector h(const Vector &m, const Matrix &D) const override;
    Matrix G(const Vector &m, const Matrix &D) const override;
    Matrix F(const Vector &m, const Matrix &D) const override;


protected:
    Matrix dadx(const Vector &x) const override;
    Matrix dadv(const Vector &x) const override;
    Matrix dbdx(const Vector &x) const override;
    Matrix dbdw(const Vector &x) const override;

    void loadParams() override;

protected:
    double m_omega, m_alpha, m_beta;

};


} // end Tasks::Discrete

} // end Tasks

#endif // DISCRETE_VANDERPOL_LINEAR_H
