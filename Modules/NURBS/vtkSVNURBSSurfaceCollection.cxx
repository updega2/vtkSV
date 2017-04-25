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

#include "vtkSVNURBSSurfaceCollection.h"
#include "vtkObjectFactory.h"

#include "vtkSVNURBSSurface.h"

vtkStandardNewMacro(vtkSVNURBSSurfaceCollection);

/**
 * Add a NURBS object to the list.
 */
void vtkSVNURBSSurfaceCollection::AddItem(vtkSVNURBSSurface *ds)
{
  this->vtkCollection::AddItem(ds);
}

/**
 * Get the next NURBS object in the list.
 */
vtkSVNURBSSurface *vtkSVNURBSSurfaceCollection::GetNextItem()
{
  return static_cast<vtkSVNURBSSurface *>(this->GetNextItemAsObject());
}

/**
 * Get the ith NURBS object in the list.
 */
vtkSVNURBSSurface *vtkSVNURBSSurfaceCollection::GetItem(int i)
{
  return static_cast<vtkSVNURBSSurface *>(this->GetItemAsObject(i));
}

/**
 * Reentrant safe way to get an object in a collection. Just pass the
 * same cookie back and forth.
 */
vtkSVNURBSSurface *vtkSVNURBSSurfaceCollection::GetNextDataObject(vtkCollectionSimpleIterator &cookie)
{
  return static_cast<vtkSVNURBSSurface *>(this->GetNextItemAsObject(cookie));
}
