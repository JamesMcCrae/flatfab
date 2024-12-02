flatfab
=======

Build
-----

### Mac OSX

- Download the Eigen library (3.2.2 stable): http://eigen.tuxfamily.org/index.php and the tar.gz and place the "Eigen" directory into the flatfab source directory
- "Release" build using Qt Creator to create flatfab.app
- May have to add the $QTDIR/bin directory to paths, do: sudo nano /etc/paths and add the line /Users/YOURNAMEHERE/Qt.5.3.2/5.3/clang_64/bin then do a control+O and enter to save then control+X to exit 
- Use icon file flatfab_icon.png, load it up and press command+A then command+C, then for the .app file in the Release directory do a control+click and "get info", then click the icon and press command+V to paste the graphic
- Terminal command: macdeployqt flatfab.app -dmg
- Repeat the icon setting process above for the flatfab.dmg file - this is the file to distribute


### Linux

You need Qt5 (tested on Ubuntu 24.10 with QtCreator 13; any previous releae should work as well)

- install `libeigen3-dev` and `libglu1-mesa-dev`
- open `FlatFab.pro` from QtCreator
- press Build and enjoy!
