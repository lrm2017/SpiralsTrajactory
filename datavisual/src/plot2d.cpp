#include "plot2d.hpp"
#include <iostream>

QSplineSeries* Plot2D::plotxy(Eigen::MatrixXd)
{
    QSplineSeries* line1 = new QSplineSeries();
    for(double x=0; x<10; x+=0.1)
    {
       line1->append(x,cos(x));
    }

//    line1->append(a,b);
    line1->setName("sin(x)");
    return line1;
}

QString Plot2D::_parseLegendPos(QString cmd)
{
    if (!cmd.size())
        return QString();

    QRegExp equalsRegex("(\\=)"); //RegEx for '='
    QStringList wordList = cmd.split(equalsRegex);

    // step1: trim spaces
    for (int i = 0; i < wordList.size(); i++)
        wordList[i] = wordList[i].trimmed();

    // step 2: delete empty strings
    QMutableStringListIterator it(wordList); // pass list as argument
    while (it.hasNext())
    {
        if (!it.next().size())
            it.remove();
    }

    int i = -1;
    QStringList::iterator iter = std::find(wordList.begin(), wordList.end(), "loc");
    if (iter != wordList.end())
        i = iter - wordList.begin(); // if "loc" is in the list, i has the index

    // ok, "label" exists && there's another word after it on the list
    if (i >= 0 && i+1 < wordList.size())
        return wordList[i+1];

    return QString();
}

void Plot2D::_parseLegend()
{
    if( !_legend.size())
        return;
    QRegExp equalsRegex("(\\=)"); //RegEx for '='   正则表达式
    QStringList wordList = _legend.split(equalsRegex);

    // step1: trim spaces
    for (int i = 0; i < wordList.size(); i++)
        wordList[i] = wordList[i].trimmed();

    // step 2: delete empty strings
    QMutableStringListIterator it(wordList); // pass list as argument
    while (it.hasNext())
    {
        if (!it.next().size())
            it.remove();
    }

    int i = -1;
    QStringList::iterator iter = std::find(wordList.begin(), wordList.end(), "label");
    if (iter != wordList.end())
        i = iter - wordList.begin(); // if "loc" is in the list, i has the index

    // ok, "label" exists && there's another word after it on the list
    if (i >= 0 && i+1 < wordList.size())
        _legend = wordList[i+1];
    else
        _legend.clear();
}

void Plot2D::show()
{
#if (DEBUG > 0) && (DEBUG < 2)
        qDebug() << "show(): " << _title;
#endif
        if (!_seriesVec.size())
        {
            qCritical() << "show()!!! Must set the data with plot() before show().";
            return;
        }

        /* Customize chart title */

        QFont font;
        font.setPixelSize(12);
        font.setWeight(QFont::Bold);
        _chart->setTitleFont(font);
        _chart->setTitle(_title);

        //TODO: investigate detaching the legend for custom positioning
        //https://doc.qt.io/qt-5/qtcharts-legend-example.html

        if (_legendPos.size())
        {
            if (_legendPos == "lower center")      // 9
                _chart->legend()->setAlignment(Qt::AlignBottom);
            else if (_legendPos == "upper center") // 8
                _chart->legend()->setAlignment(Qt::AlignTop);
            else if (_legendPos == "center right") // 7
                _chart->legend()->setAlignment(Qt::AlignRight);
            else if (_legendPos == "center left") // 6
                _chart->legend()->setAlignment(Qt::AlignLeft);
            else
            {
                qCritical() << "show()!!!" << _legendPos << " is not a valid legend position.";
                _chart->legend()->setAlignment(Qt::AlignBottom);
            }
        }

        if (_legend.size())
            _chart->legend()->setVisible(true);
        else
            _chart->legend()->setVisible(false);

        /* Customize X, Y axis and categories */

#if (DEBUG > 1) && (DEBUG < 3)
        qDebug() << "show(): xrange [" << _xMin << "," << _xMax << "] " <<
                     " yrange [" << _yMin << "," << _yMax << "]";
#endif

        QPen axisPen(Qt::black); // default axis line color and width
        axisPen.setWidth(1);

        QtCharts::QValueAxis* axisX = NULL;
        QtCharts::QCategoryAxis* categoryX = NULL;
        if (_showXticks == SHOW_TICK || _showXticks == HIDE_TICK)
        {
            bool add = true;
            if (_xAxisBottom){
                axisX = dynamic_cast<QtCharts::QValueAxis*>(_xAxisBottom);
                add = false;
            }else{
                axisX = new QtCharts::QValueAxis;
            }
            axisX->setGridLineVisible(_enableGrid);
            axisX->setTitleText(_xLabel);
            axisX->setLinePen(axisPen);
            axisX->setRange(_xMin, _xMax);
            axisX->setTickCount(_xTickCount);
            if (!_customLimits) axisX->applyNiceNumbers();
            if(add)
                _chart->addAxis(axisX, Qt::AlignBottom);
            _xAxisBottom = axisX;
            if (_showXticks == HIDE_TICK)
                axisX->setLabelsVisible(false);
        }
        else if (_showXticks == SHOW_CUSTOM_TICK)
        {
            categoryX = new QtCharts::QCategoryAxis();
            categoryX->setGridLineVisible(_enableGrid);
            categoryX->setLinePen(axisPen);

            if (_showXticks && _xTicks.size())
                for (int i = 0; i < _xTicks.size(); i++)
                {
#if (DEBUG > 1) && (DEBUG < 3)
                    qDebug() << "show(): xtick[" << i << "]=(" << _xTicks[i].second << " , " <<
                                _xTicks[i].first << ")";
#endif
                    categoryX->append(_xTicks[i].first, _xTicks[i].second);
                }

            categoryX->setRange(_xMin, _xMax);
            categoryX->setTickCount(_xTicks.size());
            if(_xAxisBottom)
                _chart->removeAxis(_xAxisBottom);
            _chart->addAxis(categoryX, Qt::AlignBottom);
            _xAxisBottom = categoryX;
        }

        QtCharts::QValueAxis* axisY = NULL;
        QtCharts::QCategoryAxis* categoryY = NULL;
        if (_showYticks == SHOW_TICK || _showYticks == HIDE_TICK) // this is the default ticks setup
        {
            bool add = true;
            if(_yAxisLeft){
                axisY = dynamic_cast<QtCharts::QValueAxis*>(_yAxisLeft);
                add = false;
            }else{
                axisY = new QtCharts::QValueAxis;
            }

            axisY->setGridLineVisible(_enableGrid);
            axisY->setTitleText(_yLabel);
            axisY->setLinePen(axisPen);
            axisY->setRange(_yMin, _yMax);
            axisY->setTickCount(_yTickCount);
            if (!_customLimits) axisY->applyNiceNumbers();
            if(add)
                _chart->addAxis(axisY, Qt::AlignLeft);
            _yAxisLeft = axisY;

            if (_showYticks == HIDE_TICK)
                axisY->setLabelsVisible(false);
        }
        else if (_showYticks == SHOW_CUSTOM_TICK) // this is for user defined ticks
        {
            categoryY = new QtCharts::QCategoryAxis();
            categoryY->setGridLineVisible(_enableGrid);
            categoryY->setLinePen(axisPen);

            if (_showYticks && _yTicks.size() > 0)
                for (int i = 0; i < _yTicks.size(); i++)
                {
#if (DEBUG > 1) && (DEBUG < 3)
                    qDebug() << "show(): ytick[" << i << "]=(" << _yTicks[i].second << " , " <<
                                _yTicks[i].first << ")";
#endif
                    categoryY->append(_yTicks[i].first, _yTicks[i].second);
                }

            categoryY->setRange(_yMin, _yMax);
            categoryY->setTickCount(_yTicks.size());
            if(_yAxisLeft)
                _chart->removeAxis(_yAxisLeft);
            _chart->addAxis(categoryY, Qt::AlignLeft);
            _yAxisLeft = categoryY;
        }

        /* Other possible customizations such as margins and background color */
        // Remove (fat) exterior margins from QChart
        _chart->layout()->setContentsMargins(0, 0, 0, 0);
        _chart->setBackgroundRoundness(0);

        /* Add series of data */
        _chart->removeAllSeries();
        for(auto i : _seriesVec)
        {
            auto itr = i.second;
            _chart->addSeries(itr.get());

            if (_showXticks == SHOW_TICK || _showXticks == HIDE_TICK)
            {
                itr->attachAxis(axisX);
            }
            else if (_showXticks == SHOW_CUSTOM_TICK)
            {
                itr->attachAxis(categoryX);
            }

            if (_showYticks == SHOW_TICK || _showYticks == HIDE_TICK)
            {
                itr->attachAxis(axisY);
            }
            else if (_showYticks == SHOW_CUSTOM_TICK)
            {
                itr->attachAxis(categoryY);
            }
        }

        _chartView->setRenderHint(QPainter::Antialiasing);
        _chartView->resize(600, 600);

        // Take a screenshot of the widget before it's destroyed
        // so it can be saved later, when savefig() is invoked after show().

        _chartView->show();

        // This loop blocks execution & waits for the window to be destroyed.
        // However, is this chart is supposed to be a real widget, then do none of this.
        if (!_isWidget)
        {
            _pixmap = _chartView->grab();
            _chartView->setAttribute(Qt::WA_DeleteOnClose); // This deletes _chartView!
            QEventLoop loop;    //
            QObject::connect(_chartView, SIGNAL(destroyed()), &loop, SLOT(quit()));
            loop.exec();

//            // _chartView was automatically deleted. Release all the other resources!
//            //_seriesVec.clear();
        }

#if (DEBUG > 0) && (DEBUG < 2)
        qDebug() << "show(): -----" ;
#endif

//    if( !_seriesVec.size())
//    {
//        qCritical() << "show()!!! Must set the data with plot() before show().";
//        return;
//    }

//    int sz = _seriesVec.size();
//    qInfo("number:%d",sz);
//    QFont font;
//    font.setPixelSize(12);
//    font.setWeight(QFont::Bold);
//    _chart->setTitleFont(font);
//    _chart->setTitle(_title);

//    if( _legend.size())
//        _chart->legend()->setVisible(true);
//    else
//        _chart->legend()->setVisible(false);

//    _chart->removeAllSeries();
//    int cnt = 0;
//    for(auto i:_seriesVec)
//    {
//        auto itr=i.second;
//        _chart->addSeries(itr.get());
//        qInfo("line:%d use_cnt:%d",cnt++,itr.use_count());
//    }

////    for(auto i:_serVec)
////    {
////        auto itr=i.second;
////        _chart->addSeries(itr);
////        qInfo("line:%d",cnt++);
////    }

//    _chartView->setRenderHint(QPainter::Antialiasing);
//    _chartView->resize(600,400);
//    _chartView->show();



//    if (!_isWidget)
//    {

//        _pixmap = _chartView->grab();
//        _chartView->setAttribute(Qt::WA_DeleteOnClose); // This deletes _chartView!

//        QEventLoop loop;
//        QObject::connect(_chartView, SIGNAL(destroyed()), &loop, SLOT(quit()));
//        loop.exec();

//        // _chartView was automatically deleted. Release all the other resources!
//        //_seriesVec.clear();
//    }
}

void Plot2D::axis(QString cmd)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "axis(): cmd=" << cmd;
#endif
    if (cmd == "off")
    {
        _showYticks = _showXticks = HIDE_TICK;
    }
    else if (cmd == "xoff")
    {
        _showXticks = HIDE_TICK;
    }
    else if (cmd == "yoff")
    {
        _showYticks = HIDE_TICK;
    }
    else
    {
        qCritical() << "axis()!!! options are 'off', 'xoff' and 'yoff'.";
        return;
    }
}

/* axis(): gets the current axes limits [xMin, xMax, yMin, yMax].
 */
void Plot2D::axis(qreal* xMin, qreal* xMax, qreal* yMin, qreal* yMax)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "axis(): _xMin=" << _xMin << " _xMax=" << _xMax << " _yMin=" << _yMin << " _yMax=" << _yMax;
#endif
    *xMin = _xMin;
    *xMax = _xMax;
    *yMin = _yMin;
    *yMax = _yMax;
}
/* axis(): sets the viewport of the axis by a list of [xMin, xMax, yMin, yMax].
 */
void Plot2D::axis(const qreal& xMin, qreal xMax, const qreal& yMin, const qreal& yMax) //设置坐标的x,y范围，输入常数
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "axis(): xMin=" << xMin << " xMax=" << xMax << " yMin=" << yMin << " yMax=" << yMax;
#endif
    _xMin = xMin;
    _xMax = xMax;
    _yMin = yMin;
    _yMax = yMax;
    _customLimits = true;
}
/* xlim(): sets the x limits of the current axes.
 */
void Plot2D::xlim(const qreal& xMin, const qreal& xMax)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "xlim(): xMin=" << xMin << " xMax=" << xMax;
#endif
    _xMin = xMin;
    _xMax = xMax;
    _customLimits = true;
}

/* ylim(): sets the x limits of the current axes.
 */
void Plot2D::ylim(const qreal& yMin, const qreal& yMax)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "ylim(): yMin=" << yMin << " yMax=" << yMax;
#endif
    _yMin = yMin;
    _yMax = yMax;
    _customLimits = true;
}

/* title(): defines the title of the chart.
 */
void Plot2D::title(QString string)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "title(): string=" << string;
#endif
    _title = string;
}

/* xlabel(): defines the label displayed below the x axis.
 */
void Plot2D::xlabel(QString label)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "xlabel(): label=" << label;
#endif
    _xLabel = label;
}

/* ylabel(): defines the label displayed to the left of the y axis.
 */
void Plot2D::ylabel(QString label)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "ylabel(): label=" << label;
#endif
    _yLabel = label;
}
/* legend(): defines the position of the legend label inside the chart.
 */
void Plot2D::legend()
{
    _legendPos = "lower center";
}

/* legend(): defines the position of the legend label inside the chart.
 */
void Plot2D::legend(QString cmd)
{
    _legendPos = _parseLegendPos(cmd);
}

/* grid(): enables or disables the background grid of the chart.
 */
void Plot2D::grid(bool status)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "grid(): status=" << status;
#endif
    _enableGrid = status;
}
/* savefig(): saves the chart displayed by show() as an image on the disk.
 */
void Plot2D::savefig(QString filename)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "savefig(): filename=" << filename;
#endif
    if (!_pixmap.isNull())
        _pixmap.save(filename);
    else if(_isWidget && _chartView){
        _pixmap = _chartView->grab();
    }
}
/* xticks(): sets the x-limits of the current tick locations and labels.
 */
void Plot2D::xticks(const Eigen::ArrayXf& values, const QVector<QString>& labels)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "xticks(): values.sz=" << values.rows() << " labels.sz=" << labels.size();
#endif
    if (values.rows() == 0 && labels.size() == 0)
    {
        _showXticks = HIDE_TICK;
        return;
    }

    if (values.rows() != labels.size())
    {
        qCritical() << "xticks(): the amount of values and labels must match!";
        return;
    }

    for (int i = 0; i < values.rows(); i++)
        _xTicks.push_back(QPair<QString, qreal>(labels[i], values[i]));

    _showXticks = SHOW_CUSTOM_TICK;

#if (DEBUG > 1) && (DEBUG < 3)
    qDebug() << "xticks(): xticks.sz=" << _xTicks.size();
    for (int i = 0; i < _xTicks.size(); i++)
        qDebug() << "\t" << _xTicks[i].second << " = " << _xTicks[i].first;
#endif
}

/* yticks(): sets the y-limits of the current tick locations and labels.
 */
void Plot2D::yticks(const Eigen::ArrayXf& values, const QVector<QString>& labels)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "yticks(): values.sz=" << values.rows() << " labels.sz=" << labels.size();
#endif
    if (values.rows() == 0 && labels.size() == 0)
    {
        _showYticks = HIDE_TICK;
        return;
    }

    if (values.rows() != labels.size())
    {
        qCritical() << "PlotLib::xticks(): the amount of values and labels must match!" ;
        return;
    }

    for (int i = 0; i < values.rows(); i++)
        _yTicks.push_back(QPair<QString, qreal>(labels[i], values[i]));

    _showYticks = SHOW_CUSTOM_TICK;

#if (DEBUG > 1) && (DEBUG < 3)
    qDebug() << "yticks(): yticks.sz=" << _yTicks.size();
    for (int i = 0; i < _yTicks.size(); i++)
        qDebug() << "\t" << _yTicks[i].second << " = " << _yTicks[i].first;
#endif
}

/* locator_params(): reduce or increase the amount of ticks for each axis.
 */
void Plot2D::locator_params(QString axis, int nbins)
{
#if (DEBUG > 0) && (DEBUG < 2)
    qDebug() << "locator_params(): axis=" << axis << " nbins=" << nbins;
#endif
    if (axis == "x")
    {
        _xTickCount = nbins;
    }
    else if (axis == "y")
    {
        _yTickCount = nbins;
    }
    else if (axis == "both")
    {
        _yTickCount = _xTickCount = nbins;
    }
    else
    {
        qCritical() << "locator_params(): '" << axis << "' is not a valid option.";
        return;
    }
}

//-----------鼠标移动 begin------------------
void Plot2D::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton)
    {
        m_lastPoint = event->pos();
        m_isPress = true;
    }
}

void Plot2D::mouseMoveEvent(QMouseEvent *event)
{
    if (!m_coordItem)
    {
        m_coordItem = new QGraphicsSimpleTextItem(this->chart());
        m_coordItem->setZValue(5);
        m_coordItem->setPos(100, 60);
        m_coordItem->show();
    }
    const QPoint curPos = event->pos();
    QPointF curVal = this->chart()->mapToValue(QPointF(curPos));
    QString coordStr = QString("X = %1, Y = %2").arg(curVal.x()).arg(curVal.y());
    m_coordItem->setText(coordStr);

    if (m_isPress)
    {
        QPoint offset = curPos - m_lastPoint;
        m_lastPoint = curPos;
        if (!m_alreadySaveRange)
        {
            this->saveAxisRange();
            m_alreadySaveRange = true;
        }
        this->chart()->scroll(-offset.x(), offset.y());
    }
}

void Plot2D::mouseReleaseEvent(QMouseEvent *event)
{
    m_isPress = false;
    if (event->button() == Qt::RightButton)
    {
        if (m_alreadySaveRange)
        {
            this->chart()->axisX()->setRange(m_xMin, m_xMax);
            this->chart()->axisY()->setRange(m_yMin, m_yMax);
        }
    }
}

//保存原始位置
void Plot2D::saveAxisRange()
{
    QValueAxis *axisX = dynamic_cast<QValueAxis*>(this->chart()->axisX());
    m_xMin = axisX->min();
    m_xMax = axisX->max();
    QValueAxis *axisY = dynamic_cast<QValueAxis*>(this->chart()->axisY());
    m_yMin = axisY->min();
    m_yMax = axisY->max();
}
//-----------鼠标移动 end -------------------

//-----------鼠标缩放 start -----------------
void Plot2D::wheelEvent(QWheelEvent *event)
{
   const QPoint curPos = event->pos();
   QPointF curVal = this->chart()->mapToValue(QPointF(curPos));

   if (!m_alreadySaveRange)
   {
    this->saveAxisRange();
    m_alreadySaveRange = true;
   }
   const double factor = 1.5;//缩放比例
   if (m_ctrlPress)
   {//Y轴
    QValueAxis *axisY = dynamic_cast<QValueAxis*>(this->chart()->axisY());
    const double yMin = axisY->min();
    const double yMax = axisY->max();
    const double yCentral = curVal.y();

    double bottomOffset;
    double topOffset;
    if (event->delta() > 0)
    {//放大
        bottomOffset = 1.0 / factor * (yCentral - yMin);
        topOffset = 1.0 / factor * (yMax - yCentral);
    }
    else
    {//缩小
        bottomOffset = 1.0 * factor * (yCentral - yMin);
        topOffset = 1.0 * factor * (yMax - yCentral);
    }

    this->chart()->axisY()->setRange(yCentral - bottomOffset, yCentral + topOffset);
   }
   else
   {//X轴
    QValueAxis *axisX = dynamic_cast<QValueAxis*>(this->chart()->axisX());
    const double xMin = axisX->min();
    const double xMax = axisX->max();
    const double xCentral = curVal.x();

    double leftOffset;
    double rightOffset;
    if (event->delta() > 0)
    {//放大
        leftOffset = 1.0 / factor * (xCentral - xMin);
        rightOffset = 1.0 / factor * (xMax - xCentral);
    }
    else
    {//缩小
        leftOffset = 1.0 * factor * (xCentral - xMin);
        rightOffset = 1.0 * factor * (xMax - xCentral);
    }
    this->chart()->axisX()->setRange(xCentral - leftOffset, xCentral + rightOffset);
   }
}
//-----------鼠标缩放 end -------------------

//-----------监听 start---------------------
//监听事件类型
//bool Plot2D::event(QEvent *event)
//{
//    switch(event->type())
//    {
////    case QEvent::TouchBegin:
////        //accepting touch begin allows us to get touch updates
////        return true;
////        break;
//    case QEvent::Gesture:               //如果是手势事件就交给手势事件处理
//        return gestureEvent(static_cast<QGestureEvent*>(event));
//        break;
//    default:                            //默认为鼠标事件
//        break;
//    }
//    return QWidget::event(event);
//}

////处理手势事件
//bool Plot2D::gestureEvent(QGestureEvent *event)
//{
//    if (!m_alreadySaveRange)            //执行手势事件之前先保存原始位置
//    {
//        this->saveAxisRange();
//        m_alreadySaveRange = true;
//    }

//    if (QGesture *pinch = event->gesture(Qt::PinchGesture))         //如果是捏合手势，进行缩放处理
//    {
//        pinchTriggered(static_cast<QPinchGesture *>(pinch));
//        event->accept();
//    }
//    return true;
//}

////处理缩放事件
//void Plot2D::pinchTriggered(QPinchGesture *gesture)
//{
//    QPinchGesture::ChangeFlags changeFlags = gesture->changeFlags();
//    if (changeFlags & QPinchGesture::ScaleFactorChanged)
//    {
//        currentStepScaleFactor = gesture->scaleFactor();            //合计放大系数
//    }
//    if (gesture->state() == Qt::GestureFinished)
//    {
//        scaleFactor = 1;
//        scaleFactor *= currentStepScaleFactor;
//        currentStepScaleFactor = 1;
//    }

//    if(scaleFactor >= 1)
//    {
//        this->chart()->zoom(1.05);
//    }
//    else if(scaleFactor < 1)
//    {
//        this->chart()->zoom(0.95);
//    }
//    update();
//}

//-----------监听 end-----------------------
