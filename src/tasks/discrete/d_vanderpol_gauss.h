#ifndef DISCRETE_VANDERPOL_GAUSS_H
#define DISCRETE_VANDERPOL_GAUSS_H

#include "d_vanderpol_linear.h"

namespace Tasks
{

namespace Discrete
{


/*!
 * \brief Осциллятор Ван-дер-Поля (дискретный, гауссовский) для фильтров оптимальной структуры.
 */

class VanDerPolGauss : public VanDerPolLinear
{

public:
    //! \brief Конструктор.
    VanDerPolGauss();


protected:
    Vector tau(const Vector &m, const Matrix &D) const override;
    Matrix Theta(const Vector &m, const Matrix &D) const override;
    Matrix A(const Vector &m, const Matrix &D) const override;
    Matrix G(const Vector &m, const Matrix &D) const override;
    Matrix F(const Vector &m, const Matrix &D) const override;

};


} // end Tasks::Discrete

} // end Tasks


#endif // DISCRETE_VANDERPOL_GAUSS_H

