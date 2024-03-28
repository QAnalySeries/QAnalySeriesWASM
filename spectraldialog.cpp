#include "spectraldialog.h"
#include "ui_spectraldialog.h"
#include "qcustomplot.h"
#include "myMath.h"
#include "spectral.h"
#include "globals.h"

QVector<double>  specxx1, specyy1, specxx2, specyy2;
QVector<double> freq1, amp1, freq2, amp2;
int specnn, specnn2;


SpectralDialog::SpectralDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SpectralDialog)
{
//    qDebug() << "in SpectralDialog init";
    ui->setupUi(this);

    connect(ui->comboBox_SpecWindow, SIGNAL(currentIndexChanged(int)), this, SLOT(on_comboBox_Input_currentIndexChanged(int)));
    connect(ui->evolCheckbox, SIGNAL(stateChanged(int)), this, SLOT(toggleEvolCheckbox_stateChanged(int)));
    connect(ui->evolWinSize, SIGNAL(currentTextChanged(Qtring)), this, SLOT(on_evolWinSize_textChanged(QString)));
    connect(ui->evolNbWindows, SIGNAL(currentTextChanged(QString)), this, SLOT(on_evolWinSize_textChanged(QString)));
    connect(ui->horizontalSlider, SIGNAL(valueChanged(int)), this, SLOT(on_horizontalSlider_valueChanged(int)));
    connect(ui->widthProductOrPercent, SIGNAL(currentTextChanged(QString)), this, SLOT(on_widthProductOrPercent_textChanged(QString)));
    connect(ui->nbWindowsOrLags, SIGNAL(currentTextChanged(QString)), this, SLOT(on_nbWindowsOrLags_textChanged(QString)));
    connect(ui->confLevel, SIGNAL(currentTextChanged(QString)), this, SLOT(on_confLevel_textChanged(QString)));
    connect(ui->logXScaleCheckbox, SIGNAL(toggled(bool)), this, SLOT(logXScaleCheckbox_toggled(bool)));
    connect(ui->logYScaleCheckbox, SIGNAL(toggled(bool)), this, SLOT(logYScaleCheckbox_toggled(bool)));

    connect(ui->applyButton, SIGNAL(clicked()), this, SLOT(applyButton_clicked()));
    connect(ui->closeButton, SIGNAL(clicked()), this, SLOT(closeButton_clicked()));

    connect(ui->saveSpectralResults, SIGNAL(clicked()), this, SLOT(saveSpectralResults_clicked()));

    spectralDataGraph = new QCustomPlot();
    spectralAnalysisGraph = new QCustomPlot();

    connect(ui->spectralDataGraph, &QCustomPlot::mouseWheel, this, &SpectralDialog::mouseWheel);
    connect(ui->spectralAnalysisGraph, &QCustomPlot::mouseWheel, this, &SpectralDialog::mouseWheel);

    ui->spectralDataGraph->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->spectralDataGraph->addGraph(0);
    ui->spectralDataGraph->addGraph();
    ui->spectralDataGraph->xAxis->setLabel("Depth or Time");
    ui->spectralDataGraph->yAxis->setLabel("Value");

    ui->spectralAnalysisGraph->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->spectralAnalysisGraph->addGraph(0);
    ui->spectralAnalysisGraph->addGraph();
    ui->spectralAnalysisGraph->xAxis->setLabel("Frequency");
    ui->spectralAnalysisGraph->yAxis->setLabel("Power");

    ui->EvolWidget->setVisible(ui->evolCheckbox->isChecked());
    ui->evolStatText->setVisible(ui->evolCheckbox->isChecked());

    ui->tabWidget->setCurrentIndex(0);

    overSamplFactor = 4.0;

    bool old = ui->comboBox_Input->blockSignals(true);
    ui->comboBox_Input->clear();

    // initial data setup
    if (x_d.length() > 0) {
        ui->comboBox_Input->addItem("Data");
        specxx1 = x_d;
        specyy1 = y_d;
        specxx2.clear();
        specyy2.clear();
        specnn = n_d;
        specnn2 = 0;
    }

    if (x_d1.length() > 0) {
        ui->comboBox_Input->addItem("Data (aged)");
        specxx1 = x_d1;
        specyy1 = y_d;
        specxx2.clear();
        specyy2.clear();
        specnn = n_d;
        specnn2 = 0;
    }

    if (x_t.length() > 0) {
        ui->comboBox_Input->addItem("Target");
        specxx1 = x_t;
        specyy1 = y_t;
        specxx2.clear();
        specyy2.clear();
        specnn = n_t;
        specnn2 = 0;
    }

    if((x_d1.length() > 0) && (x_t.length() > 0)) {
        ui->comboBox_Input->addItem("Cross Data(aged), Target");
        specxx1 = x_d1;
        specyy1 = y_d;
        specxx2 = x_t;
        specyy2 = y_t;
        specnn = n_d;
        specnn2 = n_t;
    }

    if (ui->comboBox_Input->count() == 0) {
        QMessageBox *msgBox = new QMessageBox(this);
        msgBox->setText("Please load Data and or Target first");
        msgBox->setModal(true);
        msgBox->show();
        return;
    }

    ui->comboBox_Input->blockSignals(old);

    ui->comboBox_Input->setCurrentIndex(-1); // sets specnn, specnn2 etc., also calls initData
    ui->comboBox_Input->setCurrentIndex(0); // sets specnn, specnn2 etc., also calls initData
    return;
};

SpectralDialog::~SpectralDialog()
{
    delete ui;
}

void SpectralDialog::on_comboBox_Input_currentIndexChanged(int newIndex) {
    QString selectedText;
    selectedText = ui->comboBox_Input->itemText(newIndex);

    if (selectedText == "Data") {
        specxx1 = x_d;
        specyy1 = y_d;
        specxx2.clear();
        specyy2.clear();
        specnn = n_d;
        specnn2 = 0;
    } else if (selectedText == "Data (aged)") {
        specxx1 = x_d1;
        specyy1 = y_d;
        specxx2.clear();
        specyy2.clear();
        specnn = n_d;
        specnn2 = 0;
    } else if (selectedText == "Target") {
        specxx1 = x_t;
        specyy1 = y_t;
        specxx2.clear();
        specyy2.clear();
        specnn = n_t;
        specnn2 = 0;
    } else if (selectedText == "Cross Data(aged), Target") {
        specxx1 = x_d1;
        specyy1 = y_d;
        specxx2 = x_t;
        specyy2 = y_t;
        specnn = n_d;
        specnn2 = n_t;
    } else {
        return;
    }

    initData();
    setupDataPlot();
    ui->spectralDataGraph->replot();
    ui->spectralAnalysisGraph->replot();
}


void SpectralDialog::setupDataPlot()
{
    qDebug() << "SpectralDialog: setupDataPlot";
    ui->spectralDataGraph->graph(0)->data().clear();
    ui->spectralDataGraph->graph(0)->setData(specxx1,specyy1);
    ui->spectralDataGraph->graph(0)->rescaleAxes();
    ui->spectralDataGraph->replot();
}

void SpectralDialog::setupSpectralPlot()
{
    qDebug() << "SpectralDialog: setupSpectralPlot";
    ui->spectralAnalysisGraph->graph(0)->data().clear();
    ui->spectralAnalysisGraph->graph(0)->setPen(QPen(Qt::red));
    ui->spectralAnalysisGraph->graph(0)->setSelectable(QCP::stSingleData);
    ui->spectralAnalysisGraph->graph(0)->setData(freq1,amp1);
    ui->spectralAnalysisGraph->graph(0)->rescaleAxes();
    connect(ui->spectralAnalysisGraph, &QCustomPlot::mouseMove, this, &SpectralDialog::showPointToolTip);
    ui->spectralAnalysisGraph->replot();
}

void SpectralDialog::setupDialogScale(/* STableItem1D*	item,*/ MyMath::double_regular_array* x, CommandT command )			//	initialization is here
{
    in_x = x;

    size_t n = in_x->length();
    size_t wLength = (n > 50 ? MyMath::long_round(50+70*(log10((double)n)-1.7)) : n);	//	nb of points in each window
    size_t wNumb = (2 * n)/wLength;														//	nb of windows

    double wSize = (wLength-1) * in_x->step();											//	real size of window
    double w_over_percent = (wLength * wNumb - n)/((double)n);
    ui->evolWinSize->setText(QString::number(wSize,'f',1));
    ui->evolNbWindows->setText(QString::number(wNumb));
    QString message;
    QTextStream os(&message);
    os << "(" << wLength << " points by window, " << MyMath::long_round(100 * w_over_percent) << "% oversampling)";
    ui->evolStatText->setText(message);


    int	nb = 0;//inputItem->TableSize();//	nb of series
    if (specnn > 0){
        nb += 1;
    }
    if (specnn2 > 0){
        nb += 1;
    }

    ui->crossSpectrumCheckbox->setChecked(nb > 1 );

    bool	crossAnalysis = (command == cmd_BTukey && nb>1);
    ui->crossSpectrumCheckbox->setEnabled(crossAnalysis);
    HideAndShowCrossPhase( crossAnalysis );

    HideAndShowCommand( command );
    HideAndShowEvol();
    SetParamScale( command );
    AdjustParamNumberFromPercent();

    AdjustParamNumberFromCursor();
    AdjustParamPercentFromNumber();
    AdjustStatistics();
    SetDefaultScale();
};


double SpectralDialog::linear( double x, double x1, double x2, double y1, double y2 )	{	//	linear interpolation
    return	y1 + (y2-y1) * (x-x1) / (x2-x1);
}

double SpectralDialog::paramscale( double v, const double p0[], const double p1[] ) {					//	from [0, 1] to [0, 1]
    for (int i=1; i<=paramscale_length-1; i++) {
        if ( (p0[i-1] <= v) && (v < p0[i]) ) {
            return	linear( v, p0[i-1], p0[i], p1[i-1], p1[i] );
        }
    }
    return	0;
}

void SpectralDialog::HideAndShowEvol() {
    if (ui->evolCheckbox->checkState()) {
        ui->EvolWidget->setVisible(true);
        ui->evolStatText->setVisible(true);
    } else {
        ui->EvolWidget->setVisible(false);
        ui->evolStatText->setVisible(false);
    }
};

void SpectralDialog::HideAndShowResample()
{	if (ui->resampleFreqCheckbox->checkState()) {
        ui->freqRangeWidget->show();
    } else {
        ui->freqRangeWidget->hide();
    }
}

void SpectralDialog::SetParamScale( CommandT command ) {
    switch (command) {
    case cmd_BTukey:
        current_paramscale = BT_paramscale;
        break;
    case cmd_MaxEntropy:
        current_paramscale = ME_paramscale;
        break;
    case cmd_MTM:
        current_paramscale = MTM_paramscale;
        linked_paramscale  = false;
        break;
    }
};


MyMath::spectralwindowfunc*	SpectralDialog::GetSpectralWindow() {
    int windowTypeIndex = ui->comboBox_SpecWindow->currentIndex();

    PaneIDT windowType;
    switch (windowTypeIndex) {
        case 1:	windowType = k_BartlettWindow;	break;
        case 2:	windowType = k_TukeyWindow; break;
        case 3:	windowType = k_ParzenWindow; break;
        case 4:	windowType = k_WelchWindow;	break;
        case 0:
        default: windowType = k_squareWindow;	break;
    }

    MyMath::spectralwindowfunc*	window = NULL;
    switch (windowType) {
        case k_BartlettWindow:	window = new MyMath::Bartlett_spectralwindowfunc();		break;
        case k_TukeyWindow:		window = new MyMath::Tukey_spectralwindowfunc();		break;
        case k_ParzenWindow:	window = new MyMath::Parzen_spectralwindowfunc();		break;
        case k_WelchWindow:		window = new MyMath::Welch_spectralwindowfunc();		break;
        case k_squareWindow:
        default:				window = new MyMath::square_spectralwindowfunc();		break;
   }
   return	window;
}

long SpectralDialog::series_length()	{		//	between 2 and in_x->length()
   bool	isOk  = true;
   double	wSize = ui->evolWinSize->text().toDouble(&isOk);
   if (ui->evolCheckbox->isChecked() && isOk) {
        size_t wLength = 1 + MyMath::long_round( wSize/in_x->step() );
        return std::min( in_x->length(), std::max( (size_t)2, wLength ));
    } else {
        return in_x->length();
    }
};


void SpectralDialog::AdjustParamPercentFromNumber() {
    bool isOk = true;
    int lag = ui->nbWindowsOrLags->text().toInt(&isOk);
    double	pc  = std::max( 0.0, std::min( 1.0, lag/(double)( linked_paramscale ? series_length() : 150 ) ));
    ui->widthProductOrPercent->setText(QString::number(100*MyMath::Arondi( pc, MyMath::round_nearest, 2 )));
};


void	SpectralDialog::AdjustParamNumberFromCursor() {
    int	sv = ui->horizontalSlider->value();
    double	v  = sv/100.0;
    AdjustParam( v );
};

void	SpectralDialog::AdjustParamNumberFromPercent() {
    bool isOk = true;
    double	pc = ui->widthProductOrPercent->text().toDouble(&isOk) * series_length() / 100.0;
    if (isOk) {
          int	lag = std::max((long)2, std::min( series_length(), (long)MyMath::long_round( pc ) ));
        ui->nbWindowsOrLags->setText(QString::number(lag));
    }
};

void SpectralDialog::AdjustParamCursorFromNumber() {
    bool isOk = true;
    int lag = ui->nbWindowsOrLags->text().toInt(&isOk);
    if (isOk) {
      double	v  = lag/(double)( linked_paramscale ? series_length() : 150 );
      double	gc = paramscale( v, current_paramscale, Graph_paramscale );
      long	n = MyMath::long_round(100 * gc);
      ui->horizontalSlider->setValue(n);
    }
};

void SpectralDialog::AdjustParam( double v ) {			//	v between 0 and 1, with 0.25; 0.5; 0.75 special turning points (predefined levels)
    if (current_paramscale) {
        double	pc = paramscale( v, Graph_paramscale, current_paramscale ) * ( linked_paramscale ? series_length() : 150 );					//	v * pc = simple percent value
        long	lag = std::max((long)2, std::min( series_length(), MyMath::long_round( pc ) ));
        ui->nbWindowsOrLags->setText(QString::number(lag));
    }
};

void SpectralDialog::HideAndShowParam( bool sh ) {
//    if (sh) {
//        //GetSlider( kSpectral_Param_Slider )->Show();
//        //GetStaticText( kSpectral_Confidence_StatTxt )->Show();
//        //GetStaticText( kSpectral_Resolution_StatTxt )->Show();
//        //GetStaticText( kSpectral_Percent_StatTxt )->Show();
//        //GetStaticText( kSpectral_NbLag_StatTxt )->Show();
//        //GetEditText( kSpectral_ParamPercent_ETxt )->Show();
//        //GetEditText( kSpectral_ParamNumber_ETxt )->Show();
//    } else {
//        //GetSlider( kSpectral_Param_Slider )->Hide();
//        //GetStaticText( kSpectral_Confidence_StatTxt )->Hide();
//        //GetStaticText( kSpectral_Resolution_StatTxt )->Hide();
//        //GetStaticText( kSpectral_Percent_StatTxt )->Hide();
//        //GetStaticText( kSpectral_NbLag_StatTxt )->Hide();
//        //GetEditText( kSpectral_ParamPercent_ETxt )->Hide();
//        //GetEditText( kSpectral_ParamNumber_ETxt )->Hide();
//    }
 };

void SpectralDialog::HideAndShowConfLevel( bool sh ) {
//    if (sh) {
//        //GetEditText( kSpectral_ConfLevel_ETxt )->Show();
//        //GetStaticText( kSpectral_ConfLevel_StatTxt )->Show();
//        //GetStaticText( kSpectral_ConfLevel_StatTxt2 )->Show();
//        //GetStaticText( kSpectral_StatComment_StatTxt )->Show();
//    } else {
//        //GetEditText( kSpectral_ConfLevel_ETxt )->Hide();
//        //GetStaticText( kSpectral_ConfLevel_StatTxt )->Hide();
//        //GetStaticText( kSpectral_ConfLevel_StatTxt2 )->Hide();
//        //GetStaticText( kSpectral_StatComment_StatTxt )->Hide();
//    }
};

void SpectralDialog::HideAndShowInfSupErrors( bool sh ) {
    ui->infSupErrorBarCheckbox->setVisible(sh);
};

void SpectralDialog::HideAndShowMTMAmplitude( bool sh ) {
//    	if (sh)		//GetCheckBox(kSpectral_MTM_Amplitude_ChkBox)->Show();
//        else		//GetCheckBox(kSpectral_MTM_Amplitude_ChkBox)->Hide();
};
void SpectralDialog::HideAndShowCrossPhase( bool sh ) {
//    	if (sh)		GetCheckBox(kSpectral_CrossPhase_ChkBox)->Show();
//        else		GetCheckBox(kSpectral_CrossPhase_ChkBox)->Hide();
};

void SpectralDialog::HideAndShowWindow( bool sh ) {
    ui->windowTypeWidget->setVisible(sh);
};

void SpectralDialog::HideAndShowCommand( CommandT inCommand ) {
    switch (inCommand) {
    case cmd_Periodogram:
        break;
    case cmd_BTukey:
        break;
    case cmd_MTM:
        ui->widthProductOrPercentLabel->setText("width.ndata product");
        ui->nbWindowsOrLagsLabel->setText("nb of windows");
        break;
    case cmd_MaxEntropy:
    default:
        ui->resampleFreqCheckbox->setChecked(true);
        ui->resampleFreqCheckbox->setEnabled(false);
        break;
    }

    HideAndShowParam( 		 inCommand != cmd_Periodogram );
    HideAndShowWindow( 		 (inCommand == cmd_BTukey || inCommand == cmd_Periodogram) );
    HideAndShowConfLevel( 	 inCommand == cmd_BTukey );
    HideAndShowInfSupErrors( (inCommand == cmd_BTukey || inCommand == cmd_MTM) );
    HideAndShowMTMAmplitude( inCommand == cmd_MTM );
    HideAndShowResample();
};

void SpectralDialog::AdjustStatistics() {
    bool isOk = true;
    double level = ui->confLevel->text().toDouble(&isOk) / 100.0;
    if (!isOk) {
        return;
    }

    MyMath::spectralwindowfunc*	window = GetSpectralWindow();

    int lag = ui->nbWindowsOrLags->text().toDouble(&isOk);
    if (!isOk) {
        return;
    }

    int	n = series_length();
    double	dx = in_x->step();
    QString message;
    QTextStream s(&message);

    if (lag >= 2)
    {	double	m_over_n = lag/(double)n;
        double	bw = window->BandWidth( lag )/dx;
        window->SetBTStats( level, m_over_n );
        s << "The bandwidth is " << bw;

        if (ui->crossSpectrumCheckbox->isChecked()) {
            s << "               Non-zero coherence is higher than " << window->NonZeroCoherence();
        }
        s << "\rThe error estimation on the power spectrum is " << window->LowError()
            << " < ∆Power / Power < " << window->HighError();
    }
    delete window;
};

void SpectralDialog::SetDefaultScale() {
    bool	slicing = ui->evolCheckbox->isChecked();
    double	fc 	= 1/(2*in_x->step());
    double	stp = MyMath::Arondi( 2*fc/series_length()/overSamplFactor, MyMath::round_nearest, -1 );
    double	fm  = MyMath::Arondi( ( slicing ? fc : 0.4 * fc ), MyMath::round_nearest, 2 );

    ui->freqScaleFrom->setText(QString::number(0.0, 'g', 3));
    ui->freqScaleTo->setText(QString::number(fm, 'g', 3));
    ui->freqScaleStep->setText(QString::number(stp, 'g', 3));
}

bool SpectralDialog::GetFreqScale( double& from, double& to, double& step ) {
    bool	isOk = true;
    from = ui->freqScaleFrom->text().toDouble(&isOk);
    to = ui->freqScaleTo->text().toDouble(&isOk);
    step = ui->freqScaleStep->text().toDouble(&isOk);
    return ( isOk && (step > 0) && (to > from) );
};

double SpectralDialog::OverFactor() {
    bool	isOk = true;
    double	fc 	= 1/(2*in_x->step());
    double	stp	= ui->freqScaleStep->text().toDouble(&isOk);
    int	lag = ui->nbWindowsOrLags->text().toInt(&isOk);
    int	n_zz = MyMath::nextPower2( lag-1 ) + 1;
    return 	isOk ? overSamplFactor * fc/stp/n_zz : 1.0;
}

void SpectralDialog::toggleEvolCheckbox_stateChanged(int state) {
//    qDebug() << "toggleEvolCheckbox_stateChanged";
    HideAndShowEvol();
    SetDefaultScale();
}

void SpectralDialog::on_evolWinSize_textChanged(QString newString) {
//    qDebug() << "on_evolWinSize_textChanged";
    bool isOk = true;
    double wSize = ui->evolWinSize->text().toDouble(&isOk);
    size_t wNumb = ui->evolNbWindows->text().toInt(&isOk);
    size_t wLength = 1 + MyMath::long_round(wSize/in_x->step());
    double w_over_percent = (wLength * wNumb - in_x->length())/((double)(in_x->length()));
    QString message;
    QTextStream os = QTextStream(&message);
    os << "(" << wLength << " points by window, " << MyMath::long_round(100 * w_over_percent) << "% oversampling)";
    ui->evolStatText->setText(message);
    SetDefaultScale();
}

void SpectralDialog::on_horizontalSlider_valueChanged(int newValue) {
//    qDebug() << "on_horizontalSlider_valueChanged" << newValue;
    AdjustParamNumberFromCursor();
    AdjustParamPercentFromNumber();
    AdjustStatistics();
}

void SpectralDialog::on_widthProductOrPercent_textChanged(QString newString) {
//    qDebug() << "on_widthProductOrPercent_textChanged " << newString;
    bool oldState = ui->nbWindowsOrLags->blockSignals(true);
    if (linked_paramscale)	AdjustParamNumberFromPercent();
    if (linked_paramscale)	AdjustParamCursorFromNumber();
    ui->nbWindowsOrLags->blockSignals(oldState);
    AdjustStatistics();
};

void SpectralDialog::on_nbWindowsOrLags_textChanged(QString newString) {
//    qDebug() << "on_nbWindowsOrLags_textChanged " << newString;
    bool oldState = ui->widthProductOrPercent->blockSignals(true);
    AdjustParamCursorFromNumber();
    if (linked_paramscale)	AdjustParamPercentFromNumber();
    ui->widthProductOrPercent->blockSignals(oldState);
    AdjustStatistics();
};

void SpectralDialog::on_comboBox_SpecWindow_stateChanged(int state) {
//    qDebug() << "on_nbWindowsOrLags_textChanged " << state;
    AdjustStatistics();
};

void SpectralDialog::on_confLevel_textChanged(QString newString) {
//    qDebug() << "on_confLevel_textChanged " << newString;
    AdjustStatistics();
};

void SpectralDialog::on_resampleFreqCheckbox_toggled(bool state) {
//    qDebug() << "on_resampleFreqCheckbox_stateChanged " << state;
    HideAndShowResample();
}

void SpectralDialog::initData() {
    n = 1;
    m  = specxx1.length();
    double_vector axisX = double_vector(specxx1.length());
    double_vector axisY = double_vector(specyy1.length());

    using It = decltype(std::begin(axisX));

    QVector<double>::iterator ix = specxx1.begin();
    for (It iterx = axisX.firstPtr(); ix != specxx1.end(); ix++, iterx++) {
        *iterx = *ix;
    }

    QVector<double>::iterator iy = specyy1.begin();
    for (It itery = axisY.firstPtr(); iy != specyy1.end(); iy++, itery++) {
        *itery = *iy;
    }

    val1Darray<double_vector1*,1>	in_x(1);
    in_x(1) = &axisX;
    x = new MyMath::common_regular_scale(in_x);
    va = new double_matrix1( n, m );

    va->set_all_to(0.0);

    QString msgText;
    QTextStream out(&msgText);
//    qDebug() << "Building a new (common) scale from " << (*x)(1) << " to " << (*x)(x->length()) << " with a step of " << x->step() << " (" << x->length() << " points )";
    out << "Building a new (common) scale from " << (*x)(1) << " to " << (*x)(x->length()) << " with a step of " << x->step() << " (" << x->length() << " points )";

    QMessageBox *msgBox = new QMessageBox(this);
    msgBox->setText(msgText);
    msgBox->setModal(true);
    msgBox->show();

    MyMath::double_extra_midpoints_array	mPoints( *x );
    for (size_t i=1; i<=n; i++)
    {	MyMath::lin_interpolfunc	xy_interp( axisX, axisY);
        for (size_t k=1; k<=m; k++)
            (*va)(i,k) = xy_interp.MeanBetween( mPoints[k-1], mPoints[k] );
    }

    QString commandText;
    commandText = ui->comboBox_SpecMethod->currentText();
    if (commandText == "Periodogram") {
        inCommand = cmd_Periodogram;
    } else if (commandText == "Blackman-Tukey") {
        inCommand = cmd_BTukey;
    } else if (commandText == "Maximum Entropy") {
        inCommand = cmd_MaxEntropy;
    } else if (commandText == "Multitaper Method (MTM)") {
        inCommand = cmd_MTM;
    }
    setupDialogScale(x, inCommand);
}

void SpectralDialog::logXScaleCheckbox_toggled(bool newState) {
//    qDebug() << "on_logXScaleCheckbox_toggled " << newState;
    ui->spectralAnalysisGraph->xAxis->setScaleType(newState ? QCPAxis::stLogarithmic :QCPAxis::stLinear);
    ui->spectralAnalysisGraph->replot();
}

void SpectralDialog::logYScaleCheckbox_toggled(bool newState) {
//    qDebug() << "on_logYScaleCheckbox_toggled " << newState;
    ui->spectralAnalysisGraph->yAxis->setScaleType(newState ? QCPAxis::stLogarithmic :QCPAxis::stLinear);
    ui->spectralAnalysisGraph->replot();
}

void SpectralDialog::applyButton_clicked() {
    qDebug() << "on_applyButton_clicked";

    size_t	ser_len = series_length();
    Boolean	unselect = true;

    //		the options (put the dialog readings into the "SpectralAnalyzer")
    MyMath::SpectralAnalyzer::method spectralMethod = MyMath::SpectralAnalyzer::BTukey; //default
    switch ( inCommand ) {
      case cmd_Periodogram:	spectralMethod = MyMath::SpectralAnalyzer::Periodogram;	break;
      case cmd_BTukey:		spectralMethod = MyMath::SpectralAnalyzer::BTukey;		break;
      case cmd_MaxEntropy:	spectralMethod = MyMath::SpectralAnalyzer::MaxEntropy;	break;
      case cmd_MTM:			spectralMethod = MyMath::SpectralAnalyzer::MTM;			break;
    };

    MyMath::SpectralAnalyzer	spec( spectralMethod, x->step(), ser_len );
    spec.setWindow( GetSpectralWindow() );
    bool isOk = true;
    spec.setIParam( ui->nbWindowsOrLags->text().toInt(&isOk) );
    if (!isOk) {
        return;
    }

    spec.setRParam( ui->widthProductOrPercent->text().toDouble(&isOk) );					//	for MTM, the bandwidth.
    double	fromFreq, toFreq, stepFreq;


    GetFreqScale( fromFreq, toFreq, stepFreq );
    spec.setFScale( fromFreq, toFreq, stepFreq );
    spec.setResampleFlag( ui->resampleFreqCheckbox->isChecked() );
    spec.setOverSampling(OverFactor() );
    spec.setConfLevel( ui->confLevel->text().toDouble()/100.0 );

    spec.set_inf_sup_Output( ui->infSupErrorBarCheckbox->isChecked() );
    spec.set_MTM_Ftest_Output( inCommand == cmd_MTM );												//	always for MTM
    spec.set_MTM_Ampl_Output( false /*theDialog->GetCheckBoxValue( kSpectral_MTM_Amplitude_ChkBox )*/ );
    spec.set_CrossPhase_Output( false /*theDialog->GetCheckBoxValue( kSpectral_CrossPhase_ChkBox )*/ );

    Boolean	crossAnalysis 	= (inCommand == cmd_BTukey && ui->crossSpectrumCheckbox->isChecked() );
    Boolean	multipleOut 	= spec.inf_sup_output || spec.MTM_Ftest_output || spec.MTM_Ampl_output || crossAnalysis || (n>1);

    spec.SetUpMethod();
    if (crossAnalysis)	spec.SetUpCrossAnalysis( n );		//	=> a (n x n) symetric matrix output

    if (ui->evolCheckbox->isChecked())
    {	//	2D output
        bool isOk = true;
        double	wSize = ui->evolWinSize->text().toDouble(&isOk);
        if (!isOk) {return;}
        int 	window_size = 1 + MyMath::long_round(wSize/x->step());
        size_t 	window_nb = ui->evolNbWindows->text().toInt(&isOk);
        if (!isOk) {return;}

        double	window_step = (window_nb > 1 ? (m - window_size)*x->step()/((double)(window_nb - 1)) : 0);

        for (size_t i=1; i<=n; i++)
        {
            double_matrix*	sp = new double_matrix1( window_nb, spec.spectralMethod->x_scale.length() );
            double_vector	y( *va, double_vector::fixed_x, i-1 );
            MyMath::lin_interpolfunc	lin_y( *x, y );

            for (size_t k=1; k<=window_nb; k++)
            {	double	start_i = (*x)(1) + (k-1)*window_step;
                MyMath::double_regular_array	x0( start_i, start_i+wSize, x->step() );
                double_vector	y0( x0, lin_y );
                spec.ComputePowerSpectrum( y0 );
                sp->set_colum( k, *(spec.spectrum) );
            }

            double_vector*	out_y = new double_vector1( spec.spectralMethod->x_scale );
            double_vector*	out_x = new MyMath::double_regular_array( (*x)(1) + 0.5*wSize, (size_t)window_nb, window_step );

            qDebug() << "SpectralDialog: setupPlot";

            using It = decltype(std::begin(*out_x));
            It ifreq1 = out_x->firstPtr();
            freq1 = QVector<double>(*ifreq1);

            It iAmp1 = out_y->firstPtr();
            amp1 = QVector<double>(*iAmp1);
            if (unselect){
                unselect = false;
            }
        }

    } else { //	1D output
        for (size_t i=1; i<=n; i++)
        {
            //		the computation
            double_vector	y( *va, double_vector::fixed_x, i-1 );
            spec.ComputePowerSpectrum( y );

            //		the output
            double_vector*	out_x = new double_vector1( spec.spectralMethod->x_scale );		//&(spec.spectralMethod->x_scale);
            double_vector*	out_y = new double_vector1( spec.spectrum->firstPtr(), spec.spectrum->length());

            QVector<double> qoutx = QVector<double>(out_x->length());
            QVector<double> qouty = QVector<double>(out_y->length());

            using It = decltype(std::begin(*out_x));
            QVector<double>::iterator ix = qoutx.begin();
            for (It iterx = out_x->firstPtr(); ix != qoutx.end(); ix++, iterx++) {
                *ix = *iterx;
            }

            using It = decltype(std::begin(*out_y));
            QVector<double>::iterator iy = qouty.begin();
            for (It itery = out_y->firstPtr(); iy != qouty.end(); iy++, itery++) {
                *iy = *itery;
            }

            qDebug() << "SpectralDialog: setupAnalysisPlot";
            freq1 = QVector<double>(qoutx);
            amp1  = QVector<double>(qouty);

            delete out_x;
            delete out_y;

            ui->tabWidget->setCurrentIndex(1);

            setupSpectralPlot();

            if (unselect) {
                unselect = false;
            }

            if (spec.inf_sup_output)
            {	double_vector1*	inf_x = new double_vector1( *out_x );
                double_vector1*	inf_y = spec.power_inf;
                //LStr255 iCname = Mname;	iCname += " lower-spectrum of ";	iCname += theItem->item_1D(i)->Name();
                //STableItemInfo	info( "\ppower", iCname, spec.spectralMethod->comment().c_str(), Mname, theItem->item_1D(i) );
                //STableItem1D*	powerInfItem = new STableItem1D( inf_x, inf_y, theItem->x_Axis->Name(), info );
                //powerInfItem->y_plotScale = new logScale();
                //InsertNewItem( powerInfItem, unselect );
                if (unselect)	unselect = false;

                double_vector1*	sup_x = new double_vector1( *out_x );
                double_vector1*	sup_y = spec.power_sup;
                //LStr255 sCname = Mname;	sCname += " higher-spectrum of ";	sCname += theItem->item_1D(i)->Name();
                //STableItemInfo	info1( "\ppower", sCname, spec.spectralMethod->comment().c_str(), Mname, theItem->item_1D(i) );
                //STableItem1D*	powerSupItem = new STableItem1D( sup_x, sup_y, theItem->x_Axis->Name(), info1 );
                //powerSupItem->y_plotScale = new logScale();
                //InsertNewItem( powerSupItem, unselect );
                if (unselect)	unselect = false;
            }
            if (spec.MTM_Ftest_output)
            {
                double_vector1*	f_x = new double_vector1( *out_x );
                double_vector1*	f_y = spec.MTM_ftest;
                //LStr255 fCname = Mname;	fCname += " f-test spectrum of ";	fCname += theItem->item_1D(i)->Name();
                //STableItemInfo	info( "\psignificance", fCname, spec.spectralMethod->comment().c_str(), Mname, theItem->item_1D(i) );
                //InsertNewItem( new STableItem1D( f_x, f_y, theItem->x_Axis->Name(), info ), unselect );
                if (unselect)	unselect = false;
            }
            if (spec.MTM_Ampl_output)
            {
                double_vector1*	a_x = new double_vector1( *out_x );
                double_vector1*	a_y = spec.MTM_ampl;
                //LStr255 aCname = Mname;	aCname += " amplitude spectrum of ";	aCname += theItem->item_1D(i)->Name();
                //STableItemInfo	info( "\pamplitude", aCname, spec.spectralMethod->comment().c_str(), Mname, theItem->item_1D(i) );
                //InsertNewItem( new STableItem1D( a_x, a_y, theItem->x_Axis->Name(), info ), unselect );
                if (unselect)	unselect = false;
            }
        }

               val1Darray<double_vector*,1>&	spectra = *(spec.spectralMethod->spectra);
               if (crossAnalysis)
                   for (size_t i=1; 	i<=n; i++)
                       for (size_t j=i+1; 	j<=n; j++)
                       {
                           //		the computation
                           double_vector	y( *va, double_vector::fixed_x, i-1 );
                           double_vector	z( *va, double_vector::fixed_x, j-1 );
                           spec.ComputeCrossPowerSpectrum( y, z, *(spectra(i)), *(spectra(j)) );			//	spectra(i) = spectres "bruts"

                           //		the output

                           double_vector*	out_x = new double_vector1( spec.spectralMethod->x_scale );
                           //LStr255 Mname;	Mname += spec.spectralMethod->name();
                           //LStr255 Cname = Mname;		Cname += " coherence spectrum of ";
                           //Cname += theItem->item_1D(i)->Name();	Cname += " with ";	Cname += theItem->item_1D(j)->Name();
                           //STableItemInfo	info0( "\pcoherence", Cname, spec.spectralMethod->comment().c_str(), Mname, theItem->the_item_Table() );
                           //InsertNewItem( new STableItem1D( out_x, spec.coherence, theItem->x_Axis->Name(), info0 ), unselect );
                           if (unselect)	unselect = false;

                           if (spec.inf_sup_output)
                           {
                               double_vector1*	inf_x = new double_vector1( *out_x );
                               double_vector1*	inf_y = spec.coherence_inf;
                               //LStr255 iCname;	iCname += spec.spectralMethod->name();	iCname += " lower-coherence spectrum of ";
                               //iCname += theItem->item_1D(i)->Name();	iCname += " with ";	iCname += theItem->item_1D(j)->Name();
                               //STableItemInfo	info( "\pcoherence", iCname, spec.spectralMethod->comment().c_str(), Mname, theItem->the_item_Table() );
                               //InsertNewItem( new STableItem1D( inf_x, inf_y, theItem->x_Axis->Name(), info ), unselect );
                               if (unselect)	unselect = false;

                               double_vector1*	sup_x = new double_vector1( *out_x );
                               double_vector1*	sup_y = spec.coherence_sup;
                               //LStr255 sCname;	sCname += spec.spectralMethod->name();	sCname += " higher-coherence spectrum of ";
                               //sCname += theItem->item_1D(i)->Name();	sCname += " with ";	sCname += theItem->item_1D(j)->Name();
                               //STableItemInfo	info1( "\pcoherence", sCname, spec.spectralMethod->comment().c_str(), Mname, theItem->the_item_Table() );
                               //InsertNewItem( new STableItem1D( sup_x, sup_y, theItem->x_Axis->Name(), info1 ), unselect );
                               if (unselect)	unselect = false;
                           }

                           if (spec.CrossPhase_output)
                           {
                               double_vector1*	ph_x = new double_vector1( *out_x );
                               double_vector1*	ph_y = spec.phase;
                               //LStr255 pCname;	pCname += spec.spectralMethod->name();	pCname += " cross-phase spectrum of ";
                               //pCname += theItem->item_1D(i)->Name();	pCname += " with ";	pCname += theItem->item_1D(j)->Name();
                               MyMath::lin_interpol_y_cyclicfunc* cyclicFunc = new MyMath::lin_interpol_y_cyclicfunc( -MyMath::Pi, MyMath::Pi, *ph_x, *ph_y );
                               //STableItemInfo	info( "\pphase", pCname, spec.spectralMethod->comment().c_str(), Mname, theItem->the_item_Table() );
                               //InsertNewItem( new STableItem1D( ph_x, ph_y, cyclicFunc, theItem->x_Axis->Name(), info ), unselect );
                               if (unselect)	unselect = false;

                               if (spec.inf_sup_output)
                               {
                                   double_vector1*	inf_x = new double_vector1( *out_x );
                                   double_vector1*	inf_y = spec.phase_inf;
                                   //LStr255 iCname;	iCname += spec.spectralMethod->name();	iCname += " lower-phase spectrum of ";
                                   //iCname += theItem->item_1D(i)->Name();	iCname += " with ";	iCname += theItem->item_1D(j)->Name();
                                   //STableItemInfo	info( "\pphase", iCname, spec.spectralMethod->comment().c_str(), Mname, theItem->the_item_Table() );
                                   //InsertNewItem( new STableItem1D( inf_x, inf_y, theItem->x_Axis->Name(), info ), unselect );
                                   if (unselect)	unselect = false;

                                   double_vector1*	sup_x = new double_vector1( *out_x );
                                   double_vector1*	sup_y = spec.phase_sup;
                                   //LStr255 sCname;	sCname += spec.spectralMethod->name();	sCname += " higher-phase spectrum of ";
                                   //sCname += theItem->item_1D(i)->Name();	sCname += " with ";	sCname += theItem->item_1D(j)->Name();
                                   //STableItemInfo	info1( "\pphase", sCname, spec.spectralMethod->comment().c_str(), Mname, theItem->the_item_Table() );
                                   //InsertNewItem( new STableItem1D( sup_x, sup_y, theItem->x_Axis->Name(), info1 ), unselect );
                                   if (unselect)	unselect = false;
                               }
                           }
                       }

    }		//1D output
}

void SpectralDialog::closeButton_clicked() {
//    qDebug() << "on_closeButton_clicked";
    SpectralDialog::close();
}

void SpectralDialog::closeEvent( QCloseEvent * event) {
    qDebug() << "SpectralDialog::closeEvent";
    ui->closeButton->click();
    event->accept();
}

void SpectralDialog::saveSpectralResults_clicked()
{
    int freq1Count = freq1.count();
    int amp1Count = amp1.count();

    if (freq1Count == 0 || amp1Count == 0 || freq1Count != amp1Count ) {
        QString msgText;
        QTextStream out(&msgText);
        out << "No Frequency results yet: Perform Spectral Analysis first!";
        QMessageBox *msgBox = new QMessageBox(this);
        msgBox->setText(msgText);
        msgBox->setModal(true);
        msgBox->show();
        return;
    }

    QString outString;
    QTextStream stream(&outString);
    stream << "frequency" << "\t" << "power" << Qt::endl;
    for (int i=0; i<freq1.count(); i++) {
        stream << freq1[i] << "\t" << amp1[i] << Qt::endl;
    }
//    qDebug() << "saveSpectralResults_clicked output data";
    QFileDialog::saveFileContent(outString.toUtf8(), "spectrum.txt");
}

void SpectralDialog::showPointToolTip(QMouseEvent *event)
{

    if (freq1.count() == 0 || amp1.count() == 0) return;
    int it;

    double y = ui->spectralAnalysisGraph->yAxis->pixelToCoord(event->pos().y());
    double x = ui->spectralAnalysisGraph->xAxis->pixelToCoord(event->pos().x());
    it = ui->spectralAnalysisGraph->graph(0)->findBegin(x,true);
    ui->spectralAnalysisGraph->xAxis->setLabel(QString("Frequency. Pointer: freq: %1, period: %2, amp: %3").arg(freq1[it]).arg(1.0/freq1[it]).arg(amp1[it]));
    //setToolTip(QString("freq: %1, period: %2, amp: %3").arg(x).arg(1.0/x).arg(y));
    ui->spectralAnalysisGraph->replot();

}

void SpectralDialog::mouseWheel()
{
    QList<QCPAxis *> axlist1;
    QList<QCPAxis *> axlist2;

    if(QApplication::keyboardModifiers() & Qt::ShiftModifier) {
        axlist1 << ui->spectralAnalysisGraph->xAxis <<  ui->spectralAnalysisGraph->xAxis2 << ui->spectralAnalysisGraph->yAxis <<  ui->spectralAnalysisGraph->yAxis2 ;
        axlist1 << ui->spectralDataGraph->xAxis <<  ui->spectralDataGraph->xAxis2 << ui->spectralDataGraph->yAxis <<  ui->spectralDataGraph->yAxis2;
    } else {
        axlist1.clear();
        axlist2.clear();
        axlist1 << ui->spectralAnalysisGraph->xAxis <<  ui->spectralAnalysisGraph->xAxis2;
        axlist2 << ui->spectralDataGraph->xAxis <<  ui->spectralDataGraph->xAxis2 ;

    }

    ui->spectralAnalysisGraph->axisRect()->setRangeDragAxes(axlist1);
    ui->spectralAnalysisGraph->axisRect()->setRangeZoomAxes(axlist1);

    ui->spectralDataGraph->axisRect()->setRangeDragAxes(axlist2);
    ui->spectralDataGraph->axisRect()->setRangeZoomAxes(axlist2);

}
