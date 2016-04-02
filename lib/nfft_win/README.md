nfft_matlab
===========

64 bit Windows Matlab NFFT Binaries

About
=====

Exactly what it says on the tin. Precompiled 64 bit Windows binaries for the Matlab bindings to the NFFT Library. For more information about the NFFT library see: http://www-user.tu-chemnitz.de/~potts/nfft/

Installation
============

Clone the archive into your MATLAB folder, add the new folders to your path. Try and run nfft/simple_test.m or nfsft/simple_test_nfsft.m to see if everything is working.

Possible Issues
===============

The binaries were compiled on 64 bit Windows 7, using the mingw64-x86_64-gcc provided by cygwin. The Matlab mex header files were from version 2013a, your mileage may vary trying to run these binaries on other versions of Matlab.

Compilation
===========

If you want to try and compile them yourself, I have provided a guide at: http://jpowell.co.uk/?p=77

Good luck.

Legal Information & Credits
---------------------------
Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts

This software was written by Jens Keiner, Stefan Kunis and Daniel Potts.
It was developed at the Mathematical Institute, University of
Luebeck, and at the Faculty of Mathematics, Chemnitz University of Technology.

NFFT3 is free software. You can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. If not stated otherwise, this applies to all files contained in this
package and its sub-directories. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
