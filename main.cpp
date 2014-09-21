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

    QPalette darkPalette;
    darkPalette.setColor(QPalette::Window, QColor(53,53,53));
    darkPalette.setColor(QPalette::WindowText, Qt::white);
    darkPalette.setColor(QPalette::Base, QColor(25,25,25));
    darkPalette.setColor(QPalette::AlternateBase, QColor(53,53,53));
    darkPalette.setColor(QPalette::ToolTipBase, Qt::white);
    darkPalette.setColor(QPalette::ToolTipText, Qt::white);
    darkPalette.setColor(QPalette::Text, Qt::white);
    darkPalette.setColor(QPalette::Button, QColor(53,53,53));
    darkPalette.setColor(QPalette::ButtonText, Qt::white);
    darkPalette.setColor(QPalette::BrightText, Qt::red);
    darkPalette.setColor(QPalette::Link, QColor(0, 192, 192));

    darkPalette.setColor(QPalette::Highlight, QColor(0, 192, 192));
    darkPalette.setColor(QPalette::HighlightedText, Qt::black);


    qApp->setPalette(darkPalette);

    //qApp->setStyleSheet("QToolTip { color: #ffffff; background-color: #2a82da; border: 1px solid white; }");
    qApp->setStyleSheet("QToolButton:checked { background-color: #00c0c0; color: #fff;}"
                        "QWebView { background-color: #fff; }");

    // ----------


    MainWindow w;
    w.show();
    
    return a.exec();

    QWebSettings::clearMemoryCaches();

}
