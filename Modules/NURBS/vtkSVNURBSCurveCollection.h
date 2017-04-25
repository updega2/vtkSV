/*=========================================================================
 *
 * Copyright (c) 2014-2015 The Regents of the University of California.
 * All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *=========================================================================*/

/**
 *  \class vtkSVNURBSCurveCollection
 *  \brief This is a wrapper around vtkCollection for NURBS
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVNURBSCurveCollection_h
#define vtkSVNURBSCurveCollection_h

#include "vtkCollection.h"
#include "vtkSVNURBSModule.h"

class vtkSVNURBSCurve;

class VTKSVNURBS_EXPORT vtkSVNURBSCurveCollection : public vtkCollection
{
public:
  vtkTypeMacro(vtkSVNURBSCurveCollection,vtkCollection);

  /**
   * Add a NURBS object to the list.
   */
  void AddItem(vtkSVNURBSCurve *ds);

  /**
   * Get the next NURBS object in the list.
   */
  vtkSVNURBSCurve *GetNextItem();

  /**
   * Get the ith NURBS object in the list.
   */
  vtkSVNURBSCurve *GetItem(int i);

  /**
   * Reentrant safe way to get an object in a collection. Just pass the
   * same cookie back and forth.
   */
  vtkSVNURBSCurve *GetNextDataObject(vtkCollectionSimpleIterator &cookie);

  static vtkSVNURBSCurveCollection *New();

protected:
  vtkSVNURBSCurveCollection() {}
  ~vtkSVNURBSCurveCollection() {}


private:
  // hide the standard AddItem from the user and the compiler.
  void AddItem(vtkObject *o) { this->vtkCollection::AddItem(o); };

  vtkSVNURBSCurveCollection(const vtkSVNURBSCurveCollection&);
  void operator=(const vtkSVNURBSCurveCollection&);
};

#endif
// VTK-HeaderTest-Exclude: vtkCollection.h
