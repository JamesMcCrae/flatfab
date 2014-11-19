flatfab
=======

Mac OSX build notes:

- "Release" build using Qt Creator to create flatfab.app 
- Terminal command: macdeployqt flatfab.app -dmg
- Use icon file flatfab_icon.png, requres setting both .dmg file and contained .app file to read/write mode
- Download the Eigen library (3.2.2 stable): http://eigen.tuxfamily.org/index.php and the tar.gz and place the "Eigen" directory into the flatfab source directory
