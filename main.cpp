#include <QApplication>
#include "mainwindow.h"

int main(int argc, char *argv[])
{

    int aa_samples = 16;

    if (argc > 1) {
        for (int i=1; i<argc; ++i) {
            QString eacharg(argv[i]);
            if (QString::compare(eacharg, "-ms", Qt::CaseInsensitive) == 0 && i < argc-1) {
                QString eacharg2(argv[i+1]);
                aa_samples = eacharg2.toInt();
                ++i;
            }
        }
    }

    QApplication a(argc, argv);   

    // Setting up Multi-sampling
    QGLFormat glf = QGLFormat::defaultFormat();
    glf.setSampleBuffers(aa_samples > 0);
    glf.setSamples(aa_samples);
    QGLFormat::setDefaultFormat(glf);



    // --------- New Colours

    qApp->setStyle(QStyleFactory::create("Fusion"));


//    // Dark palette
//    QPalette darkPalette;
//    darkPalette.setColor(QPalette::Window, QColor(53,53,53));
//    darkPalette.setColor(QPalette::WindowText, Qt::white);
//    darkPalette.setColor(QPalette::Base, QColor(25,25,25));
//    darkPalette.setColor(QPalette::AlternateBase, QColor(53,53,53));
//    darkPalette.setColor(QPalette::ToolTipBase, QColor(25,25,25));
//    darkPalette.setColor(QPalette::ToolTipText, QColor(0, 192, 192));
//    darkPalette.setColor(QPalette::Text, Qt::white);
//    darkPalette.setColor(QPalette::Button, QColor(53,53,53));
//    darkPalette.setColor(QPalette::ButtonText, Qt::white);
//    darkPalette.setColor(QPalette::BrightText, Qt::red);
//    darkPalette.setColor(QPalette::Link, QColor(0, 192, 192));

//    darkPalette.setColor(QPalette::Highlight, QColor(0, 192, 192));
//    darkPalette.setColor(QPalette::HighlightedText, Qt::black);


//    qApp->setPalette(darkPalette);

//    //qApp->setStyleSheet("QToolTip { color: #ffffff; background-color: #2a82da; border: 1px solid white; }");
//    qApp->setStyleSheet("QToolButton {color: #aaa;}"
//                        "QToolButton:hover {color: #fff;}"
//                        "QToolButton:checked { background-color: #00c0c0; color: #fff; }"
//                        "QWebView { background-color: #fff; }");
////    qApp->setStyleSheet("QToolButton:checked { background-color: #191919; color: #0DD;}"
////                        "QWebView { background-color: #fff; }");


    // Light Palette


    QPalette lightPalette;
    lightPalette.setColor(QPalette::Window, QColor(250,250,250));
    lightPalette.setColor(QPalette::WindowText, QColor(100,100,100));
    lightPalette.setColor(QPalette::Base, Qt::white);
    lightPalette.setColor(QPalette::AlternateBase, QColor(136,136,136));
    lightPalette.setColor(QPalette::ToolTipBase, QColor(100,100,100));
    lightPalette.setColor(QPalette::ToolTipText, Qt::white);
    lightPalette.setColor(QPalette::Text, QColor(136,136,136));
    lightPalette.setColor(QPalette::Button, Qt::white);
    lightPalette.setColor(QPalette::ButtonText, QColor(34, 192, 36));
    lightPalette.setColor(QPalette::BrightText, Qt::green);
    lightPalette.setColor(QPalette::Link, QColor(34, 192, 36));

    lightPalette.setColor(QPalette::Highlight, QColor(39, 231, 65));
    lightPalette.setColor(QPalette::HighlightedText, Qt::black);

//    lightPalette.setColor(QPalette::Light, Qt::black);
//    lightPalette.setColor(QPalette::Midlight, Qt::black);
//    lightPalette.setColor(QPalette::Mid, Qt::black);
//    lightPalette.setColor(QPalette::Dark, Qt::black);
//    lightPalette.setColor(QPalette::Shadow, Qt::black);


    qApp->setPalette(lightPalette);

    qApp->setStyleSheet("QToolButton {color: #aaa;}"
                        "QToolButton:hover {color: #22c024;}"
                        "QToolButton:checked { background-color: #22c024; color: #fff;}"
                        "QGroupBox::title {subcontrol-origin: margin; subcontrol-position: top center; padding: 0 3px; }");

    // ----------


    MainWindow w;
    w.show();
    
    return a.exec();

    QWebSettings::clearMemoryCaches();

}
