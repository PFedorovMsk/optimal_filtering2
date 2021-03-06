#ifndef ADDITIONAL_SETTINGS_WIDGET_H
#define ADDITIONAL_SETTINGS_WIDGET_H

#include "src/math/linear_algebra.h"
#include <QGroupBox>
#include <QRadioButton>
#include <QVBoxLayout>


/*!
 \brief Класс виджета, предаставляющего дополнительные настройки для численных методов.
*/

class AdditionalSettingsWidget : public QGroupBox
{
    Q_OBJECT

public:
    //! \brief Конструктор.
    explicit AdditionalSettingsWidget(QWidget *parent = nullptr);

    //! \brief Деструктор.
    ~AdditionalSettingsWidget();


private slots:
    /*!
     \brief Задействует метод SVD для псевдообращения матриц.
     \details Слот. Реакция на сигнал toggled(bool) элемента m_radioPinvSvd.
    */
    void onPinvSvdToggled(bool checked);

    /*!
     \brief Задействует метод Гревиля для псевдообращения матриц.
     \details Слот. Реакция на сигнал toggled(bool) элемента m_radioPinvGreville.
    */
    void onPinvGrevilleToggled(bool checked);


private:
    /*!
     \brief Загружает шрифты и устанавливает их параметры (начертание, размер и т.д.).
     \see FontManager.
    */
    void loadFonts();

    //! \brief Инициализирует управляющие элементы и связывает их сигналы с нужными слотами.
    void initControls();

    //! \brief Устанавливает расположение всех элементов на виджете.
    void initLayouts();


private:
    QRadioButton *m_radioPinvSvd;
    QRadioButton *m_radioPinvGreville;
};

#endif // ADDITIONAL_SETTINGS_WIDGET_H
