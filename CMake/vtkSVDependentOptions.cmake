# Copyright (c) 2014-2015 The Regents of the University of California.
# All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if(VTKSV_BUILD_MODULE_SEGMENTATION)
  message(STATUS "Forcing Module ThirdParty VMTK to be built because Segmentation module is being built")
  set(VTKSV_BUILD_THIRDPARTY_VMTK ON CACHE BOOL "Force VMTK ON" FORCE)
  message(STATUS "Forcing Module Misc to be built because Segmentation module is being built")
  set(VTKSV_BUILD_MODULE_MISC ON CACHE BOOL "Force Misc ON" FORCE)
  message(STATUS "Forcing Module NURBS to be built because Segmentation module is being built")
  set(VTKSV_BUILD_MODULE_NURBS ON CACHE BOOL "Force NURBS ON" FORCE)
  message(STATUS "Forcing Module Parameterization to be built because Segmentation module is being built")
  set(VTKSV_BUILD_MODULE_PARAMETERIZATION ON CACHE BOOL "Force NURBS ON" FORCE)
endif()

if(VTKSV_BUILD_MODULE_PARAMETERIZATION)
  message(STATUS "Forcing Module Misc to be built because Parameterization module is being built")
  set(VTKSV_BUILD_MODULE_MISC ON CACHE BOOL "Force Misc ON" FORCE)
  message(STATUS "Forcing Module NURBS to be built because Parameterization module is being built")
  set(VTKSV_BUILD_MODULE_NURBS ON CACHE BOOL "Force NURBS ON" FORCE)
endif()

if(VTKSV_BUILD_MODULE_NURBS)
  message(STATUS "Forcing Module Misc to be built because NURBS module is being built")
  set(VTKSV_BUILD_MODULE_MISC ON CACHE BOOL "Force Misc ON" FORCE)
endif()

