#include "mainwindow.h"
#include <iostream>
#include <QApplication>

QApplication *g_app = nullptr;
MainWindow *g_appWindow = nullptr;

int main(int argc, char *argv[])
{
    QCPRange range = QCPRange(0.0, 1.0);
    QCPAxisTicker ticker = QCPAxisTicker();

    g_app = new QApplication(argc, argv);
    g_appWindow = new MainWindow();

    QFile file(":/Resources/QAnalySeries.icns");
//    qDebug() << file.exists();
    QIcon icon(":/Resources/QAnalySeries.icns");
//    qDebug() << "Icon: " << icon.availableSizes();
    g_appWindow->setWindowIcon(icon);
    g_appWindow->show();
    return g_app->exec();
}
