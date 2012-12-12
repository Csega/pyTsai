To build on windows for Python 2.5
I followed the instructions here
http://boodebr.org/main/python/build-windows-extensions

I have included:
A copy of the web page above lest it dissapears "build-windows-extensions.html"

I followed the instructions for "The MinGW way" in the html file above then
Once you have done all that you simply follow the instructions in the doc folder as you would on linux but using the windows dos promt instead.

However, python needed to be called as c:\python25\python.exe as the path was not added by the installer like this:
c:\python25\python.exe setup.py build
etc.
AND
c:\python25\Lib\distutils\version.py
raised an error which was to do with parsing mingw's version number I figured out through googling, there was a ptch on python.org to fix it but
I was not authorised to download it so I hacked it to make it work.
The modified version.py file is included USE AT YOUR OWN RISK, the part I modified is between
##START Lobo
AND
##END Lobo
I am well aware this is an ugly hack but hey it got the job done and I dont intend to compile any other python modules.
