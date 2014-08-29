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

    MainWindow w;
    w.show();
    
    return a.exec();

    QWebSettings::clearMemoryCaches();

}
