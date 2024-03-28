#ifndef SPECTRALDIALOG_H
#define SPECTRALDIALOG_H

#include <QDialog>
#include <qcustomplot.h>
#include "myMath.h"
#include "globals.h"

namespace Ui {
class SpectralDialog;
}

class SpectralDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SpectralDialog(QWidget *parent = nullptr);
    ~SpectralDialog();
    QCustomPlot *spectralDataGraph;
    QCustomPlot *spectralAnalysisGraph;

    void setupDataPlot();
    void setupSpectralPlot();


private slots:
    void mouseWheel();
    void saveSpectralResults_clicked();
    void toggleEvolCheckbox_stateChanged(int);
    void applyButton_clicked();
    void closeButton_clicked();
    void on_comboBox_Input_currentIndexChanged(int newIndex);
    void on_evolWinSize_textChanged(QString newString);
    void on_horizontalSlider_valueChanged(int);
    void on_widthProductOrPercent_textChanged(QString newString);
    void on_nbWindowsOrLags_textChanged(QString newString);
    void on_comboBox_SpecWindow_stateChanged(int);
    void on_confLevel_textChanged(QString);
    void on_resampleFreqCheckbox_toggled(bool);
    void logXScaleCheckbox_toggled(bool);
    void logYScaleCheckbox_toggled(bool);
    void showPointToolTip(QMouseEvent *event);


private:
    Ui::SpectralDialog *ui;
    void setupDialogScale(/* STableItem1D*	item,*/ MyMath::double_regular_array* x, CommandT command );
    double linear( double x, double x1, double x2, double y1, double y2 );
    double paramscale( double v, const double p0[], const double p1[] );
    void HideAndShowEvol();
    void HideAndShowResample();
    void SetParamScale( CommandT command );
    void AdjustParamPercentFromNumber();
    void AdjustParamNumberFromCursor();
    void AdjustParamNumberFromPercent();
    void AdjustParamCursorFromNumber();
    void AdjustParam( double v);
    void HideAndShowParam( bool sh );
    void HideAndShowConfLevel( bool sh );
    void HideAndShowInfSupErrors( bool sh );
    void HideAndShowMTMAmplitude( bool sh );
    void HideAndShowCrossPhase( bool sh );
    void HideAndShowWindow( bool sh );
    void HideAndShowCommand( CommandT inCommand );
    void AdjustStatistics();
    void SetDefaultScale();
    void initData();
    void closeEvent( QCloseEvent *);


//private:
    static const int	paramscale_length = 5;

    static constexpr double	Graph_paramscale[] 	= { 0, 0.25, 0.5, 0.75, 1.0 };			//	list of special predefined values
    static constexpr double	BT_paramscale[]     = { 0, 0.1, 0.3, 0.5, 1.0 };			//	list of special predefined values
    static constexpr double	ME_paramscale[] 	= { 0, 0.05, 0.1, 0.3, 1.0 };			//	list of special predefined values
    static constexpr double	MTM_paramscale[] 	= { 0.0, 0.02, 0.03, 0.06, 1.0 };		//	list of special predefined values

    double	overSamplFactor;
    const double*	current_paramscale = NULL;
    bool			linked_paramscale = true;			//	false for MTM (2 independent params, except when adjusted by cursor)
    MyMath::double_regular_array* 	in_x ;

    MyMath::common_regular_scale*	x = NULL;
    double_matrix1*					va = NULL;
    CommandT inCommand = cmd_None;
    size_t			n = 0;
    size_t			m = 0;

public:
    MyMath::spectralwindowfunc* GetSpectralWindow();
    bool GetFreqScale( double& from, double& to, double& step );
    double OverFactor();
    long	series_length();
};




#endif // SPECTRALDIALOG_H
