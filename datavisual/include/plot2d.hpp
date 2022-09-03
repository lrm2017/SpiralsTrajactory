#ifndef PLOT2D_HPP
#define PLOT2D_HPP

#include <type_traits>  //泛型编程头文件
#include <memory>
#include <iostream>

#include <QDebug>
#include <QtCharts>
#include "Eigen/Dense"
#include <QDebug>
#include <QEventLoop>
#include <QGraphicsLayout>
#include <QPair>

#include <QMouseEvent> //鼠标缩放图表
#include <QGraphicsSimpleTextItem>

#pragma once

#ifndef PLT_ARG_NAMESPACE
#define PLT_ARG_NAMESPACE
#endif

#ifdef _MSC_VER
#if _MSC_VER >= 1900
#define PLT_CONSTEXPR constexpr
#else
#define PLT_CONSTEXPR
#endif
#else
#define PLT_CONSTEXPR constexpr //验证是否为常量表达式 例如1+2+3就是常量表达式 函数的返回值被constexpr修饰过，函数只能只有一个return语句
#endif

#define MO_KEYWORD_INPUT(name, type)                                                                                                  \
namespace tag{                                                                                                                        \
    struct name {                                                                                                                     \
        typedef type        Type;                                                                                                     \
        typedef const Type& ConstRef;                                                                                                 \
        typedef Type&       Ref;                                                                                                      \
        typedef ConstRef    StorageType;                                                                                              \
        typedef const void* VoidType;                                                                                                 \
        template <typename T>                                                                                                         \
        static PLT_CONSTEXPR bool AllowedType() { return std::is_same<Type, T>::value; }                                              \
        static VoidType GetPtr(const Type& arg) {                                                                                     \
            return &arg;                                                                                                              \
        }                                                                                                                         \
        template <class T>                                                                                                            \
        static VoidType GetPtr(const T& arg) {                                                                                        \
            (void)arg;                                                                                                                \
            return nullptr;                                                                                                           \
        }                                                                                                                             \
    };                                                                                                                                \
}                                                                                                                                     \
namespace PLT_ARG_NAMESPACE {                                                                                                         \
    static kwargs::TKeyword<tag::name>& name = kwargs::TKeyword<tag::name>::instance;                                                 \
}

namespace kwargs {
    struct TaggedBase {};
    template <class Tag>
    struct TaggedArgument : public TaggedBase {
        typedef Tag TagType;
        explicit TaggedArgument(typename Tag::StorageType val)  //explicit显示赋值，避免隐式赋值
            : arg(&val) {
        }

        typename Tag::VoidType get() const {
            return arg;
        }

    protected:
        typename Tag::VoidType arg;
    };

    template <class Tag>
    struct TKeyword {
        static TKeyword     instance;
        TaggedArgument<Tag> operator=(typename Tag::StorageType data) {
            return TaggedArgument<Tag>(data);
        }
    };
    template <class T>
    TKeyword<T> TKeyword<T>::instance;
}


template <class Tag, bool Infer = false>
typename Tag::VoidType GetKeyImpl() {   //如果没有找到与def类型一致的类，返回0
    return 0;
}

template <class Tag, bool Infer = false, class T, class... Args>
typename std::enable_if<std::is_base_of<kwargs::TaggedBase, T>::value, typename Tag::VoidType>::type
    GetKeyImpl(const T& arg, const Args&... args) {
    return std::is_same<typename T::TagType, Tag>::value ? const_cast<void*>(arg.get()) : const_cast<void*>(GetKeyImpl<Tag, Infer, Args...>(args...));  //const_cast去掉const 属性
}

template <class Tag, bool Infer = false, class T, class... Args>
typename std::enable_if<!std::is_base_of<kwargs::TaggedBase, T>::value, typename Tag::VoidType>::type
    GetKeyImpl(const T& arg, const Args&... args) {
#ifdef __GNUC__
    //static_assert(CountType<typename Tag::Type>(arg, args...) <= 1, "Cannot infer type when there are multiple variadic Params with desired type");
#endif
    return Tag::template AllowedType<T>() ? // This infers the type
        Tag::GetPtr(arg)
        : const_cast<void*>(GetKeyImpl<Tag, Infer, Args...>(args...));
}

template <class Tag, bool Infer = false, class... Args>
typename Tag::ConstRef GetKeywordInputDefault(typename Tag::ConstRef def, const Args&... args) {

    const void* ptr = GetKeyImpl<Tag, Infer>(args...);
    if (ptr)
        return *static_cast<const typename Tag::Type*>(ptr);
    return def;
}



MO_KEYWORD_INPUT(x, Eigen::ArrayXf)
MO_KEYWORD_INPUT(y, Eigen::ArrayXf)
MO_KEYWORD_INPUT(marker, QString)
MO_KEYWORD_INPUT(label, QString)
MO_KEYWORD_INPUT(color, QColor)
MO_KEYWORD_INPUT(linewidth, quint32)
MO_KEYWORD_INPUT(alpha, qreal)
MO_KEYWORD_INPUT(edgecolor, QColor)
MO_KEYWORD_INPUT(markersize, qreal)

/* Debug control */

#define DEBUG               0      // level 0: debug msgs are disabled
                                   // level 1: print method calls
                                   // level 2: print method calls + data

/* Global definitions */
#define SHOW_TICK           1
#define HIDE_TICK           2
#define SHOW_CUSTOM_TICK    4

#define DEFAULT_LEGEND      ""
#define DEFAULT_MARKER      "-"
#define DEFAULT_ALPHA       1.0f
#define DEFAULT_COLOR       "none"
#define DEFAULT_EDGECOLOR   "none"
#define DEFAULT_LINEW       2
#define DEFAULT_MARKERSZ    6.0f

template<typename Derived>
struct is_array_expression     //判断是否为array矩阵表达式 is_base_of<A,B>  A是否是B的基类
    : std::is_base_of<Eigen::ArrayBase<std::decay_t<Derived> >, std::decay_t<Derived> > //decay_t<>返回退化的引用类型
{};

template<typename Derived>
struct is_matrix_expression     //判断是否为matrix矩阵表达式 is_base_of<A,B>  A是否是B的基类
    : std::is_base_of<Eigen::MatrixBase<std::decay_t<Derived> >, std::decay_t<Derived> >
{};

QT_CHARTS_USE_NAMESPACE

class Plot2D:public QChartView{


public:

    Plot2D(QChartView* chartView = NULL)
        :_chart(NULL), _chartView(chartView)
    {
#if (DEBUG > 0) && (DEBUG < 2)
        qDebug() << "Madplotlib(): isWidget=" << isWidget;
#endif

        _xAxisTop = _xAxisBottom = _yAxisLeft = _yAxisRight = nullptr;
        _chart = new QChart();
//        _chartView = new QChartView(_chart);
        if(_chartView == NULL)
            _chartView = new QChartView(_chart);
        else
        {
            _chartView->setChart(_chart);
            _isWidget = true;
        }

        _enableGrid = false;
        _customLimits = false;
        _xMin = _xMax = _yMin = _yMax = 0;

        _showXticks = _showYticks = SHOW_TICK;  //刻度间隔
        _xTickCount = 7;
        _yTickCount = 5;

        _colorIdx = 0;
        _colors = {QColor(0x1f77b4),QColor(0xff7f0e),QColor(0x2ca02c),QColor(0xd62728),QColor(0x9467bd),
                   QColor(0x8c564b),QColor(0xe377c2),QColor(0x7f7f7f),QColor(0xbcbd22),QColor(0x17becf)};
    }

    void axis(QString cmd); //隐藏坐标
    void axis(qreal* xMin, qreal* xMax, qreal* yMin, qreal* yMax);  //获取坐标的x,y范围，输入变量
    void axis(const qreal& xMin, qreal xMax, const qreal& yMin, const qreal& yMax); //设置坐标的x,y范围，输入常数
    void xlim(const qreal& xMin, const qreal& xMax);    //x轴限制
    void ylim(const qreal& yMin, const qreal& yMax);    //y轴限制
    void title(QString string);
    void xlabel(QString label);
    void ylabel(QString label);
    void legend();
    void legend(QString cmd);       //标注的位置
    void grid(bool status);         //是否显示网格
    void savefig(QString filename); //将图表保存为图片，输入文件名
    void xticks(const Eigen::ArrayXf& values, const QVector<QString>& labels);  //显示刻度以及刻度的标签
    void yticks(const Eigen::ArrayXf& values, const QVector<QString>& labels);
    void locator_params(QString axis, int nbins);

    void setMarker(QString marker){ _marker = marker;};
    void setMarkerSize(qreal mSize) { _markersize = mSize;}
    void setAlpha(qreal alpha) { _alpha = alpha;}
    void setColor(QColor color) { _color = color;}
    void setLineWidth(quint32 lw ) { _linewidth = lw;}
    void setEdgeColor(QColor egColor) { _edgecolor = egColor;}
    void setLabel(QString label) { _label = label;}

    void saveAxisRange();   //保存坐标区域，用于复位

    /* plot(): called when user needs to put data on a chart.
     * x: an array that stores x axis values.
     * y: an array that stores y axis values.
     * marker: chart types. "-" for line plot and "o" for scatter plot.
     * alpha: defines the transparency level of the color.
     * color: defines the color used to draw the data.
     * edgecolor: defines the edge color of "o" marker.
     * linewidth: defines the width of the pen used to draw "-" marker.
     * markersize: defines the size of "o" marker.
     */
    template<class T, class ... Args>
    typename std::enable_if<is_array_expression<T>::value>::type plot(const Eigen::ArrayXf& x, const T& y, const Args&... args)
    {
        plotXY(x, y, args...);
    }

    template<class ... Args>
    void plotXY(const Eigen::ArrayXf& x, const Eigen::ArrayXf& y, const Args&... args)
    {
        const QString marker  = GetKeywordInputDefault<tag::marker>(DEFAULT_MARKER, args...);
        const QString label  = GetKeywordInputDefault<tag::label>(DEFAULT_LEGEND, args...);
        auto alpha = GetKeywordInputDefault<tag::alpha>(DEFAULT_ALPHA, args...);
        const QColor color = GetKeywordInputDefault<tag::color>(DEFAULT_COLOR, args...);
        auto linewidth = GetKeywordInputDefault<tag::linewidth>(DEFAULT_LINEW, args...);
        const QColor edgecolor = GetKeywordInputDefault<tag::edgecolor>(DEFAULT_EDGECOLOR, args...);
        auto markersize = GetKeywordInputDefault<tag::markersize>(DEFAULT_MARKERSZ, args...);


#if (DEBUG > 0) && (DEBUG < 2)
        qDebug() << "plot(x,y): marker:" << marker << " alpha:" << alpha <<
                    " color:" << color << " edgecolor:" << edgecolor <<
                    " linewidth:" << linewidth << " markersize:" << markersize;
#endif

        if (marker != "-" && marker != "--" && marker != "." && marker != "o" && marker != "s")
        {
            qCritical() << "plot(x,y): unknown marker '" << marker << "'.";
            return;
        }

        if (x.rows() != y.rows())
        {
            qCritical() << "plot(x,y): x.sz=" << x.cols() << " != y.sz=" << y.rows();
            exit(-1);
        }

        if (x.rows() == 0 || y.rows() == 0)
        {
            qCritical() << "plot(x,y): axis size must be > 0 but it is x.sz=" << x.rows() << " y.sz=" << y.rows();
            exit(-1);
        }

        // Make a copy because it's show() who setup these things
        _legend = label;

        // find min and max values to define the range of the X axis
        qreal xMin = x.minCoeff();
        qreal xMax = x.maxCoeff();
        if (xMin < _xMin)
            _xMin = xMin;
        if (xMax > _xMax)
            _xMax = xMax;

        // find min and max values to stablish the range of the Y axis
        // however, if a new series brings more xtreme values, we need to respect that!
        qreal yMin = y.minCoeff();
        qreal yMax = y.maxCoeff();
        if (yMin < _yMin)
            _yMin = yMin;
        if (yMax > _yMax)
            _yMax = yMax;

#if (DEBUG > 1) && (DEBUG < 3)
        qDebug() << "plot(x,y): xrange [" << _xMin << "," << _xMax <<
                     "]  yrange [" << _yMin << "," << _yMax << "]";
#endif

        std::shared_ptr<QtCharts::QXYSeries> series;
//        auto itr = _seriesVec.find(_legend);  //想不通作者为啥要写这个
        /*if(itr != _seriesVec.end()){
//            series = itr->second;
        }else*/{
            if (marker == "o" || marker == "s") // it's a scatter plot!
            {
    #if (DEBUG > 1) && (DEBUG < 3)
                qDebug() << "plot(x,y): scatter plot";
    #endif
                QtCharts::QScatterSeries* s = new QtCharts::QScatterSeries();
                s->setMarkerSize(markersize); // symbol size

                if (marker == "o")
                    s->setMarkerShape(QtCharts::QScatterSeries::MarkerShapeCircle);

                if (marker == "s")
                    s->setMarkerShape(QtCharts::QScatterSeries::MarkerShapeRectangle);

                series.reset((QtCharts::QXYSeries*)s);
                series->setUseOpenGL(true);
            }
            else // draw line
            {
    #if (DEBUG > 1) && (DEBUG < 3)
                qDebug() << "plot(x,y): line plot";
    #endif
                series.reset(new QtCharts::QLineSeries());
                series->setUseOpenGL(true);
            }
        }
        // Call a string parser! Ex: "label=Trump Tweets" becomes "Trump Tweets"
        _parseLegend();
        if (_legend.size())
        {
#if (DEBUG > 1) && (DEBUG < 3)
            qDebug() << "plot(x,y): label=" << _legend;
#endif
            series->setName(_legend);
        }

        if(series->count() == x.rows()){
            for (int i = 0; i < x.rows(); i++)
                series->replace(i, x[i], y[i]);
        }else{
            series->clear();
            for (int i = 0; i < x.rows(); i++)
            {
#if (DEBUG > 1) && (DEBUG < 3)
                qDebug() << "plot(x,y): x[" << i << "]=" << x[i] << " y[" << i << "]=" << y[i];
#endif
                series->append(x[i], y[i]);
            }
        }

        // Customize series color and transparency
        QColor fillColor = color;
        if (fillColor == DEFAULT_COLOR)
            fillColor = _colors[_colorIdx++];
        fillColor.setAlphaF(alpha);

        QPen pen = series->pen();
        pen.setWidth(linewidth);

        if (marker == "o" || marker == "s")
        {
            if (edgecolor == DEFAULT_EDGECOLOR)
            {
                pen.setColor(fillColor);        // outline should be invisible
#if (DEBUG > 1) && (DEBUG < 3)
                qDebug() << "plot(x,y): fillColor=" << fillColor;
#endif
            }
            else
            {
#if (DEBUG > 1) && (DEBUG < 3)
                qDebug() << "plot(x,y): edgecolor=" << edgecolor;
#endif
                QColor edgeColor = edgecolor;
                edgeColor.setAlphaF(alpha);
                pen.setColor(edgeColor);       // create circle outline
            }
        }
        else if (marker == "--")
        {
            pen.setStyle(Qt::DashLine);
            pen.setColor(fillColor);
        }
        else if (marker == ".")
        {
            pen.setStyle(Qt::DotLine);
            pen.setColor(fillColor);
        }
        else // marker == "-"
        {
            pen.setColor(fillColor);
        }

        series->setPen(pen);
        series->setBrush(QBrush(fillColor));

        if (_colorIdx >= _colors.size())
            _colorIdx = 0;

        //_seriesVec.push_back(series);
        _seriesVec.insert(std::make_pair(_legend,series));


        _marker = DEFAULT_MARKER;
        _label = DEFAULT_LEGEND;
        _alpha = DEFAULT_ALPHA;
        _color = DEFAULT_COLOR;
        _linewidth = DEFAULT_LINEW;
        _edgecolor = DEFAULT_EDGECOLOR;
        _markersize = DEFAULT_MARKERSZ;

#if (DEBUG > 0) && (DEBUG < 2)
        qDebug() << "plot(x,y): -----";
#endif
    }


    QSplineSeries* plotxy(Eigen::MatrixXd);
    QSplineSeries* plotxy(Eigen::MatrixXf);

//    QSplineSeries *plot();


#if (DEBUG > 0) && (DEBUG < 2)
        qDebug() << "show(): -----" ;
#endif

    void show();    //通过调用plot()显示所有的数据
    void clear(){  _seriesVec.clear();  }   //存储曲线的列表清零

protected:
    //-----------鼠标缩放图表------------
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
//    void keyPressEvent(QKeyEvent *event);
//    void keyReleaseEvent(QKeyEvent *event);

    //以下3个为Qt Android下Qchart的缩放（单指触点时默认为鼠标点击，所以移动功能可正常使用）
//    bool event(QEvent *event) override;             //使用手势实现缩放
//    bool gestureEvent(QGestureEvent *event);
//    void pinchTriggered(QPinchGesture *gesture);
    //--------------------------------------

private:
    bool _is_marker(const QString& cmd)
    {
        if (cmd == "-" || cmd == "--" || cmd == "." || cmd == "o" || cmd == "s")
            return true;

        return false; // cmd is a label for the legend
    }

    void _check_cmds_are_good(const QString& cmd1, const QString& cmd2)
    {
        if (_is_marker(cmd1) && !_is_marker(cmd2))
            return;

        if (!_is_marker(cmd1) && _is_marker(cmd2))
            return;

        qCritical() << "_check_cmds_are_good()!!! Only one marker and one label are allowed.";
        exit(-1);
    }

    void _parseLegend();                    //legend 命令 默认是center right
    QString _parseLegendPos(QString cmd);   //legend 命令 默认是center right


    QPixmap _pixmap;                                // show() screenshots the widget, savefig() writes it on the disk
    QtCharts::QChart* _chart;                       // manages the graphical representation of the chart's series, legends & axes
    QtCharts::QChartView* _chartView;               // standalone widget that can display charts
    std::multimap<QString, std::shared_ptr<QtCharts::QXYSeries>> _seriesVec;      // every plot() creates a new series of data that is stored here
    std::multimap<QString, QXYSeries*> _serVec;  //test smart ptr different from common ptr

    bool _isWidget = false;                                // true: show() doesn't block so this can be used as widget
    QString _legend;
    QString _legendPos;

    QVector<QPair<QString, qreal> > _xTicks;       // user defined <label, endValue> ticks that replace default ticks
    QVector<QPair<QString, qreal> >  _yTicks;
    int _showXticks;                               // flag that show/hides the exhibition of ticks on the X axis
    int _showYticks;                               // flag that show/hides the exhibition of ticks on the Y axis
    int _xTickCount;                              // number of ticks displayed on the X axis
    int _yTickCount;                              // number of ticks displayed on the Y axis

    QString _title;                                 // Chart title
    QString _yLabel;                                // String that appears to the left of the Y axis
    QString _xLabel;                                // String that appears below the X axis

    QVector<QColor> _colors;                        // Vector that stores chart's predefined colours
    int _colorIdx;                                 // Every plot() increases this index so the next has a different colour

    bool _customLimits;                            // If the user has informed new limits through xlim(), ylim(), or axis()
    qreal _xMin;                                   // X axis min limit
    qreal _xMax;                                   // X axis max limit
    qreal _yMin;                                   // Y axis min limit
    qreal _yMax;                                   // Y axis max limit

    bool _enableGrid;                              // flag that show/hides the background grid
    QtCharts::QAbstractAxis* _yAxisLeft;
    QtCharts::QAbstractAxis* _yAxisRight;
    QtCharts::QAbstractAxis* _xAxisBottom;
    QtCharts::QAbstractAxis* _xAxisTop;

    QString _marker = DEFAULT_MARKER;
    QString _label = DEFAULT_LEGEND;
    qreal _alpha = DEFAULT_ALPHA;
    QColor _color = DEFAULT_COLOR;
    quint32 _linewidth = DEFAULT_LINEW;
    QColor _edgecolor = DEFAULT_EDGECOLOR;
    qreal _markersize = DEFAULT_MARKERSZ;

    //鼠标缩放图表
    QPoint m_lastPoint;
    bool m_isPress;
    bool m_ctrlPress;
    bool m_alreadySaveRange;
    double m_xMin, m_xMax, m_yMin, m_yMax;
    QGraphicsSimpleTextItem* m_coordItem;

};




#endif // PLOT2D_HPP
