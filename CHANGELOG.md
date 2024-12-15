Changelog
=========

release 0.8.4
-------------

remove all network-related code (welcome webpage and analytics),
UI refactoring, put quicktools in left toolbar

release 0.8.3
-------------

add toggle in edit menu to disable clipping with ground plane

release 0.8.2
-------------

use middle mouse button/mousewheel for more standard 3D navigation

release 0.8.1
-------------

modernize for Qt5, formatting, fix all compiler warnings

release 0.8.0
-------------

added initial support for finger joints for the surface facet generation
(teeth between adjacent faces) added ability to export slab geometry as a
2D layout (link to video) https://www.youtube.com/watch?v=b-i05v8l9nM

release 0.7.1
-------------

fixed an issue with mouse cursor position on OSX with retina display

release 0.7
-----------

stability improvement: new ear clipping algorithm performs triangulations
for sections with/without holes

release 0.6 -
-------------

more ui improvements
new interaction for specifying planes for procedural modelling operations
updated curve filtering/fitting algorithm for input strokes
potentially fixed bug on mac with cpu usage
fixed bug where "new flatfab" from menu would not first remove the
webview added filename to window title

release 0.5 -
-------------

ui improvements including
- new tool sidebar
- new transform widget

release 0.4 -
-------------

added "calibration shape" generation to make finding the right
calibration setting a cinch multisampling now a command line parameter,
the default value is 4.  e.g: flatfab.exe -ms 0 would disable
multisampling smoother animations/frame updates updated linux binary with
RPATH set in the executable (rather than relying on a script)

release 0.3
-----------

updated name
updated program icon
added antialiasing/multisampling option

release 0.2
-----------

fixed rippling effect of bezier curve fit
made 80 degrees the default rotation angle   

