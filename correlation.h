#ifndef CORRELATION_H
#define CORRELATION_H
#include <qcustomplot.h>

namespace Ui {
  class Correlation;
}

class Correlation : public QDialog
{
  Q_OBJECT

public:
  explicit Correlation(QWidget *parent = 0);
  ~Correlation();
  QCustomPlot *corrGraphic;

private:
  Ui::Correlation *ui;

};

#endif // CORRELATION_H
